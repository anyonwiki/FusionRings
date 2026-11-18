# test/properties_tests.jl
# Uses the revised test_utils.jl, test_fixtures.jl, and keyword-style runner.

using Test
using Random: seed!
using Oscar: QQBar, ZZ, order

_prop_partition(xs) = sort([sort(collect(Int, block)) for block in xs])


# ============================================================
# Comparison and witness helpers
# ============================================================

# Character rows may be reordered. Preserve multiplicity.
function _prop_check_character_rows(
  actual,
  expected,
  context;
  approximate = false,
)
  A, E = Matrix(actual), Matrix(expected)

  check_equal(
    size(A),
    size(E),
    "$context: wrong character-table shape",
  )

  size(A) == size(E) || return nothing
  remaining = collect(axes(E, 1))

  for i in axes(A, 1)
    position = findfirst(remaining) do j
      if approximate
        isapprox(A[i, :], E[j, :]; atol = 1e-8, rtol = 1e-8)
      else
        A[i, :] == E[j, :]
      end
    end

    check_true(
      position !== nothing,
      "$context: unmatched or repeated character row $i",
    )

    position === nothing && return nothing
    deleteat!(remaining, position)
  end

  return nothing
end

function _prop_check_character_identities(C, mt, context)
  A = Matrix(C)
  r = size(mt, 1)

  check_equal(
    size(A),
    (r, r),
    "$context: wrong character-table shape",
  )

  size(A) == (r, r) || return nothing

  for j in 1:r
    check_equal(
      A[j, 1],
      QQBar(1),
      "$context: row $j is not unital",
    )

    for a in 1:r, b in 1:r
      check_equal(
        A[j, a] * A[j, b],
        sum(mt[a, b, c] * A[j, c] for c in 1:r),
        "$context: character $j fails on product ($a,$b)",
      )
    end
  end

  return nothing
end

function _prop_check_dimensions(ds, mt, context)
  r = size(mt, 1)

  check_equal(
    length(ds),
    r,
    "$context: wrong dimension count",
  )

  length(ds) == r || return nothing

  check_equal(
    ds[1],
    QQBar(1),
    "$context: unit dimension is not one",
  )

  check_true(
    all(d -> isreal(d) && real(d) > 0, ds),
    "$context: dimensions are not positive real numbers",
  )

  for a in 1:r, b in 1:r
    check_equal(
      ds[a] * ds[b],
      sum(mt[a, b, c] * ds[c] for c in 1:r),
      "$context: dimension identity fails on product ($a,$b)",
    )
  end

  return nothing
end

function _prop_valid_permutation(p, r)
  return p isa AbstractVector &&
         length(p) == r &&
         sort(collect(p)) == collect(1:r) &&
         p[1] == 1
end

# Every positive equivalence result must supply a valid tensor witness.
function _prop_equivalent(a, b)
  rank(a) == rank(b) || return false

  candidates = which_permutation(a, b; force_compute = true)
  candidates === nothing && return false

  A, B = multiplication_table(a), multiplication_table(b)

  return any(
    p -> _prop_valid_permutation(p, rank(a)) && A[p, p, p] == B,
    candidates,
  )
end

function _prop_same_multiset(
  actual,
  expected;
  equal = _prop_equivalent,
)
  length(actual) == length(expected) || return false
  remaining = collect(expected)

  for item in actual
    index = findfirst(other -> equal(item, other), remaining)
    index === nothing && return false
    deleteat!(remaining, index)
  end

  return isempty(remaining)
end

function _prop_check_subring(parent, subset, child, context)
  r = rank(parent)

  valid = !isempty(subset) &&
          length(unique(subset)) == length(subset) &&
          1 in subset &&
          all(i -> 1 <= i <= r, subset)

  check_true(valid, "$context: invalid embedding indices")
  valid || return nothing

  S = sort(collect(subset))
  N = multiplication_table(parent)
  outside = setdiff(collect(1:r), S)

  check_true(
    all(iszero, N[S, S, outside]),
    "$context: product escapes the proposed subring",
  )

  check_equal_tensor(
    multiplication_table(child),
    N[S, S, S],
    "$context: wrong restricted table",
  )

  check_equal(
    labels(child),
    labels(parent)[S],
    "$context: restriction changed labels",
  )

  return nothing
end

function _prop_check_grading(parent, grades, quotient, context)
  r, q = rank(parent), rank(quotient)

  valid = length(grades) == r &&
          all(g -> 1 <= g <= q, grades)

  check_true(valid, "$context: invalid grading values")
  valid || return nothing

  check_equal(
    grades[1],
    1,
    "$context: unit has the wrong grade",
  )

  check_equal(
    sort(unique(grades)),
    collect(1:q),
    "$context: grading is not surjective",
  )

  N = multiplication_table(parent)
  Q = multiplication_table(quotient)

  for a in 1:q, b in 1:q
    check_true(
      sum(Q[a, b, :]) == 1 &&
      count(x -> x != 0, Q[a, b, :]) == 1,
      "$context: quotient is not a group fusion ring",
    )
  end

  for a in 1:r, b in 1:r, c in 1:r
    N[a, b, c] == 0 && continue

    check_equal(
      Q[grades[a], grades[b], grades[c]],
      1,
      "$context: grading fails on fusion channel ($a,$b,$c)",
    )
  end

  return nothing
end


# ============================================================
# Reference execution
# ============================================================

# Supplement spatial sampling with explicit output branches.
# Each probe adds the first matching case.
# Full mode continues to visit every case.
function _prop_reference(f::Function, filename; probes = ())
  @testset "$filename" begin
    data = load_anyonica_data(filename)
    indices = oracle_case_indices(filename, data; max_cases = 8)

    if oracle_mode() in (:representative, :spread)
      for probe in probes
        index = findfirst(
          x -> x !== nothing && probe(x),
          data["Output"],
        )

        index === nothing || push!(indices, index)
      end
    end

    for index in sort(unique(indices))
      code = data["Input"][index]
      expected = data["Output"][index]

      @testset "case $index, code $code" begin
        if expected === nothing
          @info "Unavailable reference answer" filename index code
          @test_skip false
        else
          r = ring_from_anyonica_code(code; fresh = true)
          f(r, expected)
        end
      end
    end
  end

  return nothing
end


@testset "Properties" begin

  # ==========================================================
  # Basic accessors and classification
  # ==========================================================

  maybe_testset(
    "accessors, duals, and classification";
    level = "basic",
  ) do
    for f in small_ring_fixtures()
      @testset "$(f.name)" begin
        r, e = fixture_ring(f), f.expected

        check_equal(rank(r), e.rank, "$(f.name): wrong rank")

        check_equal_tensor(
          multiplication_table(r),
          e.mt,
          "$(f.name): wrong table accessor",
        )

        check_equal(labels(r), e.labels, "$(f.name): wrong labels")

        check_equal(
          names(r)["miscellaneous"],
          [String(f.name)],
          "wrong name accessor",
        )

        check_equal(
          tex_names(r)["miscellaneous"],
          [String(f.name)],
          "wrong TeX-name accessor",
        )

        check_true(
          ismissing(anyonwiki_code(r)),
          "fresh fixture unexpectedly has a catalogue code",
        )

        dual = conjugate_element(r)

        check_equal(
          [conjugate_element(r, i) for i in 1:e.rank],
          e.duals,
          "$(f.name): wrong duals",
        )

        check_equal(
          dual.(collect(1:e.rank)),
          e.duals,
          "curried dual accessor disagrees",
        )

        check_equal(
          dual.(e.duals),
          collect(1:e.rank),
          "duality is not an involution",
        )

        expected_C = [
          Int(e.duals[i] == j)
          for i in 1:e.rank, j in 1:e.rank
        ]

        check_equal(
          conjugation_matrix(r),
          expected_C,
          "$(f.name): wrong conjugation matrix",
        )

        sd = count(i -> e.duals[i] == i, 1:e.rank)

        check_equal(
          nsdnsd(r),
          [sd, e.rank - sd],
          "wrong self-duality counts",
        )

        check_equal(nsd(r), sd, "wrong self-dual count")
        check_equal(nnsd(r), e.rank - sd, "wrong non-self-dual count")

        check_equal(
          is_commutative(r),
          e.commutative,
          "$(f.name): wrong commutativity",
        )

        check_equal(
          is_group_ring(r),
          e.group_ring,
          "$(f.name): wrong group-ring classification",
        )

        if e.group_ring
          check_equal(
            cayley_table(r),
            f.input.cayley_table,
            "wrong Cayley table",
          )
        else
          check_throws(
            () -> cayley_table(r),
            ArgumentError,
            "non-group ring accepted";
            message_contains = "group fusion ring",
          )
        end
      end
    end
  end

  maybe_testset(
    "support and multiplicity";
    level = "basic",
  ) do
    for f in small_ring_fixtures()
      @testset "$(f.name)" begin
        r = fixture_ring(f)
        N = f.expected.mt

        expected = [
          (a, b, c)
          for a in 1:rank(r)
          for b in 1:rank(r)
          for c in 1:rank(r)
          if N[a, b, c] > 0
        ]

        check_equal(
          sort(nzsc(r)),
          sort(expected),
          "wrong nonzero structure-constant coordinates",
        )

        check_equal(nnzsc(r), length(expected), "wrong support size")
        check_equal(multiplicity(r), maximum(N), "wrong multiplicity")
        check_equal(mult(r), multiplicity(r), "multiplicity alias disagrees")
      end
    end

    r = fixture_ring(fixture_rank2_m2())

    check_equal(
      multiplicity(r),
      2,
      "multiplicity-two case was lost",
    )

    check_equal(
      nnzsc(r),
      5,
      "support count incorrectly counted multiplicities",
    )
  end

  # ==========================================================
  # Dimensions and cache behavior
  # ==========================================================

  maybe_testset(
    "exact dimensions and integrality";
    level = "intermediate",
  ) do
    for f in small_ring_fixtures()
      @testset "$(f.name)" begin
        r, e = fixture_ring(f), f.expected
        ds = fpdims(r; force_compute = true)

        check_equal(
          ds,
          e.fpdims,
          "$(f.name): incorrect exact FP dimensions",
        )

        check_equal(
          fpdim(r; force_compute = true),
          e.fpdim,
          "incorrect total FP dimension",
        )

        _prop_check_dimensions(ds, e.mt, string(f.name))

        check_equal(
          sum(d*d for d in ds),
          e.fpdim,
          "sum-of-squares identity failed",
        )

        check_equal(is_integral(r), e.integral, "wrong integrality")

        check_equal(
          is_weakly_integral(r),
          e.weakly_integral,
          "wrong weak integrality",
        )

        check_approx(
          numeric_fpdims(r; force_compute = true),
          Float64.(real.(ComplexF64.(e.fpdims))),
          "numeric dimensions disagree",
        )

        check_approx(
          numeric_fpdim(r; force_compute = true),
          Float64(real(ComplexF64(e.fpdim))),
          "numeric total dimension disagrees",
        )
      end
    end
  end

  maybe_testset(
    "stored and forced dimension paths";
    level = "intermediate",
  ) do
    f = fixture_ising()
    bad = QQBar.([1, 1, 1])

    r = fusion_ring(
      copy(f.input.mt);
      frobenius_perron_dimensions = bad,
      frobenius_perron_dimension = QQBar(1),
    )

    check_equal(fpdims(r), bad, "default dimensions did not read their cache")
    check_equal(fpdim(r), QQBar(1), "default total did not read its cache")

    check_equal(
      fpdims(r; force_compute = true),
      f.expected.fpdims,
      "forced dimensions reused cache",
    )

    check_equal(
      fpdim(r; force_compute = true),
      f.expected.fpdim,
      "forced total reused nested caches",
    )

    check_approx(
      numeric_fpdims(r),
      [1.0, 1.0, 1.0],
      "numeric cache retrieval changed",
    )

    check_approx(
      numeric_fpdim(r),
      1.0,
      "numeric scalar cache retrieval changed",
    )

    check_approx(
      numeric_fpdims(r; force_compute = true),
      [1.0, sqrt(2.0), 1.0],
      "numeric force ignored",
    )

    check_approx(
      numeric_fpdim(r; force_compute = true),
      4.0,
      "numeric total force ignored",
    )

    components_only = fusion_ring(
      copy(f.input.mt);
      frobenius_perron_dimensions = bad,
    )

    check_equal(
      fpdim(components_only),
      QQBar(3),
      "default total did not use component cache",
    )

    check_equal(
      fpdim(components_only; force_compute = true),
      QQBar(4),
      "forced total reused component cache",
    )
  end

  # ==========================================================
  # Characters and formal codegrees
  # ==========================================================

  maybe_testset(
    "exact characters and formal codegrees";
    level = "intermediate",
  ) do
    for (index, f) in enumerate(small_ring_fixtures())
      @testset "$(f.name)" begin
        r, e = fixture_ring(f), f.expected

        if e.commutative
          seed!(20261118 + index)

          C = characters(r; force_compute = true)

          _prop_check_character_rows(C, e.characters, string(f.name))
          _prop_check_character_identities(C, e.mt, string(f.name))

          check_equal(
            formal_codegrees(r; force_compute = true),
            e.formal_codegrees,
            "wrong character-based codegrees",
          )

          check_equal(
            formal_codegrees(
              r;
              force_compute = true,
              usecharacters = false,
            ),
            e.formal_codegrees,
            "wrong central-element codegrees",
          )

          seed!(20261218 + index)
          numeric = numeric_characters(r; force_compute = true)

          check_true(
            numeric isa Matrix{ComplexF64},
            "numeric characters have the wrong type",
          )

          _prop_check_character_rows(
            numeric,
            ComplexF64.(e.characters),
            string(f.name);
            approximate = true,
          )
        else
          # S3 is genuinely noncommutative.
          for compute in (
            characters,
            numeric_characters,
            formal_codegrees,
          )
            check_throws(
              () -> compute(r; force_compute = true),
              ErrorException,
              "unsupported computation did not report its scope";
              message_contains = "non-commutative",
            )
          end
        end
      end
    end

    check_equal(
      numeric_characters(fixture_ring(fixture_unit())),
      ComplexF64[1.0;;],
      "rank-one numeric characters changed",
    )
  end

  maybe_testset(
    "character caches and computation controls";
    level = "intermediate",
  ) do
    f = fixture_ising()
    bad = fill(QQBar(1), 3, 3)

    r = fusion_ring(
      copy(f.input.mt);
      characters = bad,
      formal_codegrees = QQBar.([7, 7, 7]),
    )

    check_equal(
      Matrix(characters(r)),
      bad,
      "default characters ignored their cache",
    )

    check_equal(
      formal_codegrees(r),
      QQBar.([7, 7, 7]),
      "default codegrees ignored their cache",
    )

    seed!(41)

    _prop_check_character_rows(
      characters(r; force_compute = true),
      f.expected.characters,
      "forced characters",
    )

    check_equal(
      formal_codegrees(r; force_compute = true),
      f.expected.formal_codegrees,
      "forced codegrees reused stored codegrees or characters",
    )

    fresh = fixture_ring(fixture_z2())

    check_throws(
      () -> characters(fresh; force_compute = true, tries = 0),
      ErrorException,
      "exact-character controls were ignored";
      message_contains = "0 tries",
    )

    # Current numeric_characters signature takes tries and tol positionally.
    check_throws(
      () -> numeric_characters(fresh, 0, 1e-12; force_compute = true),
      ErrorException,
      "numeric attempt count was ignored";
      message_contains = "0 tries",
    )

    seed!(42)

    check_throws(
      () -> numeric_characters(fresh, 1, -Inf; force_compute = true),
      ErrorException,
      "numeric tolerance was ignored";
      message_contains = "1 tries",
    )

    seed!(43)

    exact = characters(
      fresh;
      use_numerics = false,
      force_compute = true,
    )

    _prop_check_character_rows(
      exact,
      fixture_z2().expected.characters,
      "symbolic character path",
    )
  end

  # ==========================================================
  # Subrings and embeddings
  # ==========================================================

  maybe_testset(
    "subrings, closure, and embeddings";
    level = "intermediate",
  ) do
    for f in small_ring_fixtures()
      @testset "$(f.name)" begin
        r, e = fixture_ring(f), f.expected

        subsets = sub_fusion_ring_subsets(r)

        check_equal(
          _prop_partition(subsets),
          _prop_partition(e.proper_subrings),
          "wrong subring subsets",
        )

        children = sub_fusion_rings(
          r;
          force_compute = true,
          represent_by_known = false,
        )

        check_equal(
          _prop_partition(first.(children)),
          _prop_partition(e.proper_subrings),
          "wrong materialized subrings",
        )

        for (S, child) in children
          check_true(
            is_sub_fusion_ring(r, reverse(S)),
            "subring predicate depends on index order",
          )

          _prop_check_subring(r, S, child, string(f.name))
        end

        check_equal(
          is_simple(r),
          isempty(e.proper_subrings),
          "wrong simplicity classification",
        )
      end
    end

    r = fixture_ring(fixture_z4())

    for S in (
      Int[],
      [2, 3],
      [1, 1],
      [1, 0],
      [1, 5],
      [1, 2],
    )
      check_false(
        is_sub_fusion_ring(r, S),
        "invalid subring subset $S was accepted",
      )
    end

    check_true(
      is_sub_fusion_ring(r, [3, 1]),
      "valid reordered Z4 subring was rejected",
    )

    check_false(
      is_sub_fusion_ring(fixture_ring(fixture_ising()), [1, 2]),
      "Ising sigma product escaped unnoticed",
    )

    check_equal(
      FusionRings._fusion_closure(r, [3]),
      [1, 3],
      "wrong closure of g^2",
    )

    check_equal(
      FusionRings._fusion_closure(r, [2]),
      [1, 2, 3, 4],
      "wrong closure of g",
    )

    small = fixture_ring(fixture_z2())
    injection = which_injection(small, r)

    check_true(
      injection !== nothing,
      "Z2 embedding into Z4 was missed",
    )

    if injection !== nothing
      image = [injection[i] for i in 1:2]

      check_equal(sort(image), [1, 3], "wrong Z2 image in Z4")
      check_equal(image[1], 1, "embedding moved the unit")

      check_equal_tensor(
        multiplication_table(r)[image, image, image],
        fixture_z2().expected.mt,
        "embedding does not preserve multiplication",
      )
    end

    check_equal(
      which_injection(fixture_ring(fixture_fibonacci()), r),
      nothing,
      "false Fibonacci embedding",
    )

    stale = change_properties(
      fixture_ring(fixture_v4()),
      :sub_fusion_rings => Tuple{Vector{Int},String}[],
    )

    check_true(
      isempty(sub_fusion_rings(stale)),
      "stored empty subring cache was ignored",
    )

    check_equal(
      length(sub_fusion_rings(
        stale;
        force_compute = true,
        represent_by_known = false,
      )),
      3,
      "forced subring computation reused the empty cache",
    )
  end

  # ==========================================================
  # Equivalence and automorphisms
  # ==========================================================

  maybe_testset(
    "equivalence and automorphism witnesses";
    level = "intermediate",
  ) do
    for f in small_ring_fixtures()
      @testset "$(f.name)" begin
        r, e = fixture_ring(f), f.expected
        actual = fusion_ring_automorphisms(r)

        check_equal(
          sort(actual),
          sort(e.automorphisms),
          "incomplete or duplicated automorphisms",
        )

        for p in actual
          check_true(
            _prop_valid_permutation(p, e.rank),
            "invalid automorphism permutation",
          )

          if _prop_valid_permutation(p, e.rank)
            check_equal_tensor(
              e.mt[p, p, p],
              e.mt,
              "automorphism does not preserve multiplication",
            )
          end
        end

        check_equal(
          Int(order(automorphism_group(r; force_compute = true))),
          length(e.automorphisms),
          "wrong automorphism-group order",
        )

        p = vcat(1, reverse(collect(2:e.rank)))
        relabelled = fusion_ring(e.mt[p, p, p])

        check_true(
          _prop_equivalent(r, relabelled),
          "equivalence missed a valid relabeling",
        )

        check_true(
          is_equivalent_fusion_ring(r, relabelled),
          "public equivalence predicate disagrees",
        )
      end
    end

    check_false(
      is_equivalent_fusion_ring(
        fixture_ring(fixture_z4()),
        fixture_ring(fixture_v4()),
      ),
      "Z4 and V4 were incorrectly identified",
    )
  end

  # ==========================================================
  # Adjoint rings, central series, and universal gradings
  # ==========================================================

  maybe_testset(
    "adjoint rings, central series, and universal gradings";
    level = "intermediate",
  ) do
    for f in small_ring_fixtures()
      @testset "$(f.name)" begin
        r, e = fixture_ring(f), f.expected

        S, child = adjoint_fusion_ring(
          r;
          force_compute = true,
          represent_by_known = false,
        )

        check_equal(sort(S), e.adjoint_subset, "wrong adjoint subset")
        _prop_check_subring(r, S, child, "adjoint $(f.name)")

        chain = upper_central_series(
          r;
          force_compute = true,
          represent_by_known = false,
        )

        check_equal(
          [rank(item[2]) for item in chain],
          e.upper_central_ranks,
          "wrong central-series ranks",
        )

        parent = r

        for (subset, term) in chain
          _prop_check_subring(
            parent,
            subset,
            term,
            "central-series $(f.name)",
          )
          parent = term
        end

        check_equal(
          is_nilpotent(
            r;
            force_compute = true,
            represent_by_known = false,
          ),
          e.nilpotent,
          "wrong nilpotency",
        )

        components = adjoint_irreps(
          r;
          force_compute = true,
          represent_by_known = false,
        )

        check_equal(
          _prop_partition(components),
          _prop_partition(e.grading_components),
          "wrong adjoint components",
        )

        grades, quotient = universal_grading(
          r;
          force_compute = true,
          represent_by_known = false,
        )

        _prop_check_grading(r, grades, quotient, string(f.name))

        fibers = [
          findall(==(g), grades)
          for g in 1:rank(quotient)
        ]

        check_equal(
          _prop_partition(fibers),
          _prop_partition(e.grading_components),
          "wrong universal-grading fibers",
        )
      end
    end
  end

  maybe_testset(
    "all gradings and stored grading metadata";
    level = "intermediate",
  ) do
    for (f, expected_ranks) in (
      (fixture_unit(), [1]),
      (fixture_z4(), [1, 2, 4]),
      (fixture_v4(), [1, 2, 2, 2, 4]),
      (fixture_ising(), [1, 2]),
      (fixture_s3(), [1, 2, 6]),
    )
      @testset "$(f.name)" begin
        r = fixture_ring(f)
        actual = all_gradings(r; force_compute = true)

        check_equal(
          sort([rank(item[2]) for item in actual]),
          expected_ranks,
          "wrong quotient grading ranks",
        )

        partitions = [
          _prop_partition([
            findall(==(g), grades)
            for g in 1:rank(quotient)
          ])
          for (grades, quotient) in actual
        ]

        check_equal(
          length(unique(partitions)),
          length(actual),
          "duplicate quotient gradings",
        )

        for (grades, quotient) in actual
          _prop_check_grading(
            r,
            grades,
            quotient,
            "all gradings $(f.name)",
          )
        end
      end
    end

    # UUID cache resolution intentionally uses catalogue rings.
    unit = fawc(1, 1, 0, 1)
    r = fixture_ring(fixture_z4())

    stale = change_properties(
      r,
      :all_gradings => [([1, 1, 1, 1], uuid(unit))],
    )

    check_equal(
      length(all_gradings(stale)),
      1,
      "stored grading metadata was ignored",
    )

    check_equal(
      sort([
        rank(x[2])
        for x in all_gradings(stale; force_compute = true)
      ]),
      [1, 2, 4],
      "forced gradings reused cache",
    )

    known_z2 = fawc(2, 1, 0, 1)

    cached = change_properties(
      known_z2,
      :upper_central_series => [
        (collect(1:2), uuid(known_z2)),
      ],
    )

    check_false(
      is_nilpotent(cached),
      "default nilpotency ignored stored central-series metadata",
    )

    check_equal(
      is_nilpotent(
        cached;
        force_compute = true,
        represent_by_known = false,
      ),
      true,
      "forced central-series computation reused stored metadata",
    )

    cached_grades, _ = universal_grading(
      cached;
      represent_by_known = false,
    )

    forced_grades, _ = universal_grading(
      cached;
      force_compute = true,
      represent_by_known = false,
    )

    check_equal(
      cached_grades,
      [1, 1],
      "default universal grading ignored the stored adjoint data",
    )

    check_equal(
      forced_grades,
      [1, 2],
      "forced universal grading reused the stored adjoint data",
    )
  end

  # ==========================================================
  # Commutator subsets
  # ==========================================================

  maybe_testset(
    "commutator subsets";
    level = "intermediate",
  ) do
    f = fixture_ising()
    r = fixture_ring(f)
    unit = fixture_ring(fixture_unit())

    S, child = commutator(r, ([1], unit))

    check_equal(
      S,
      [1, 3],
      "wrong Ising commutator over the unit subring",
    )

    _prop_check_subring(r, S, child, "Ising commutator")

    S2, child2 = commutator(r)

    check_equal(
      S2,
      [1, 2, 3],
      "commutator over the full ring was not full",
    )

    _prop_check_subring(r, S2, child2, "full commutator")

    check_throws(
      () -> commutator(fixture_ring(fixture_s3())),
      ErrorException,
      "noncommutative commutator did not report its scope";
      message_contains = "commutative",
    )
  end

  # ==========================================================
  # Decompositions and realization caches
  # ==========================================================

  maybe_testset(
    "decompositions and realization caches";
    level = "intermediate",
  ) do
    z2 = fixture_ring(fixture_z2())
    v4 = fixture_ring(fixture_v4())
    z4 = fixture_ring(fixture_z4())

    actual = decompositions(
      v4;
      force_compute = true,
      represent_by_known = false,
    )

    check_true(
      _prop_same_multiset(
        actual,
        [[z2, z2]];
        equal = (a, b) -> _prop_same_multiset(a, b),
      ),
      "V4 decomposition is not Z2 tensor Z2",
    )

    check_true(
      isempty(decompositions(
        z4;
        force_compute = true,
        represent_by_known = false,
      )),
      "Z4 has a false tensor decomposition",
    )

    check_throws(
      () -> decompositions(z2, "DirectSum"),
      ArgumentError,
      "unsupported decomposition kind was accepted";
      message_contains = "tensor-product",
    )

    check_true(
      FusionRings.is_equivalent_decomposition(
        [z2, z4],
        [
          fixture_ring(fixture_z4()),
          fixture_ring(fixture_z2()),
        ],
      ),
      "factor order or UUIDs changed decomposition equivalence",
    )

    check_false(
      FusionRings.is_equivalent_decomposition(
        [z2, z2],
        [z2, z4],
      ),
      "factor multiplicity was ignored",
    )

    stale = change_properties(
      v4,
      :sub_fusion_rings => Tuple{Vector{Int},String}[],
    )

    check_true(
      isempty(decompositions(
        stale;
        represent_by_known = false,
      )),
      "stored subring cache was ignored",
    )

    forced = decompositions(
      stale;
      force_compute = true,
      represent_by_known = false,
    )

    check_true(
      _prop_same_multiset(
        forced,
        [[z2, z2]];
        equal = (a, b) -> _prop_same_multiset(a, b),
      ),
      "forced decomposition reused cache",
    )

    cached_empty = change_properties(
      v4,
      :realizations => Dict{String,Any}(
        "tensor_product" => Vector{String}[],
      ),
    )

    check_true(
      isempty(decompositions(
        cached_empty;
        represent_by_known = false,
      )),
      "default decomposition ignored realization cache",
    )

    forced_realization = decompositions(
      cached_empty;
      force_compute = true,
      represent_by_known = false,
    )

    check_true(
      _prop_same_multiset(
        forced_realization,
        [[z2, z2]];
        equal = (a, b) -> _prop_same_multiset(a, b),
      ),
      "forced decomposition reused the realization cache",
    )

    for cache in (
      Dict{String,Any}(),
      Dict{String,Any}("tensor_product" => missing),
      Dict{String,Any}("other" => Any[]),
    )
      r = change_properties(
        fixture_ring(fixture_z2()),
        :realizations => cache,
      )

      check_equal(
        realizations(r)["tensor_product"],
        Vector{FusionRing}[],
        "missing realization key was not computed",
      )
    end
  end

  # ==========================================================
  # Barcode and tensor encoding
  # ==========================================================

  maybe_testset(
    "barcode and multi-digit tensor encoding";
    level = "intermediate",
  ) do
    check_equal(
      barcode(
        fixture_ring(fixture_unit());
        force_compute = true,
      ),
      ZZ(1),
      "rank-one barcode failed",
    )

    f = fixture_ising()
    p = [1, 3, 2]

    check_equal(
      barcode(fixture_ring(f); force_compute = true),
      barcode(
        fusion_ring(f.input.mt[p, p, p]);
        force_compute = true,
      ),
      "barcode depends on basis numbering",
    )

    # Valid rank-two ring x*x = 1+10x.
    mt = copy(fixture_rank2_m2().input.mt)
    mt[2, 2, 2] = 10

    expected = ZZ(0)

    for a in 1:2, b in 1:2, c in 1:2
      expected = 11*expected + mt[a, b, c]
    end

    check_equal(
      FusionRings.mult_tab_code(mt, 10),
      expected,
      "multi-digit coefficients were split into separate digits",
    )
  end

  # ==========================================================
  # Stored categorifiability metadata
  # ==========================================================

  maybe_testset(
    "categorifiability metadata retrieval";
    level = "basic",
  ) do
    r = fixture_ring(fixture_z2())

    check_true(
      categories_with_properties(r) isa AbstractDict,
      "category metadata is not a dictionary",
    )

    check_equal(
      is_categorifiable(r),
      missing,
      "fresh ring unexpectedly has a categorification verdict",
    )

    for flag in (true, false, missing)
      metadata = Dict{String,Vector{Any}}(
        "Fusion" => Any[flag, "", "test metadata"],
      )

      tagged = change_properties(
        r,
        :has_categories_with_props => metadata,
      )

      check_equal(
        is_categorifiable(tagged),
        flag,
        "stored category flag was changed",
      )
    end
  end

  # ==========================================================
  # Anyonica reference calculations
  # ==========================================================

  maybe_testset(
    "Anyonica property calculations";
    level = "reference",
  ) do
    _prop_reference(
      "fpdims.json";
      probes = (x -> any(v -> json_float(v) > 1.1, x),),
    ) do r, expected
      ds = fpdims(r; force_compute = true)

      _prop_check_dimensions(
        ds,
        multiplication_table(r),
        "reference FP dimensions",
      )

      check_numeric_vector_approx(
        ds,
        expected,
        "FP dimensions disagree with the reference",
      )
    end

    _prop_reference("fpdim.json") do r, expected
      check_numeric_approx(
        fpdim(r; force_compute = true),
        expected,
        "total FP dimension disagrees",
      )
    end

    _prop_reference(
      "is_group_ring.json";
      probes = (x -> x === true, x -> x === false),
    ) do r, expected
      check_true(expected isa Bool, "group-ring oracle is not Boolean")

      check_equal(
        is_group_ring(r),
        expected,
        "group-ring classification disagrees",
      )
    end

    _prop_reference(
      "num_nonzero_structure_constants.json",
    ) do r, expected
      check_equal(
        nnzsc(r),
        json_int(expected),
        "support size disagrees",
      )
    end

    _prop_reference(
      "sub_ring_data.json";
      probes = (isempty, x -> !isempty(x)),
    ) do r, expected
      actual = sub_fusion_rings(
        r;
        force_compute = true,
        represent_by_known = false,
      )

      check_equal(
        length(unique(_prop_partition(first.(actual)))),
        length(actual),
        "duplicate subring embeddings",
      )

      for (S, child) in actual
        _prop_check_subring(r, S, child, "reference subring")
      end

      expected_rings = [
        ring_from_anyonica_code(pair[2]; fresh = true)
        for pair in expected
      ]

      check_true(
        _prop_same_multiset(last.(actual), expected_rings),
        "subring types or multiplicities disagree",
      )
    end

    # Some exports used Anyonica's SortedRing.
    # Compare invariant quantities and validate witnesses instead of
    # assuming that both programs use identical basis numbering.

    _prop_reference(
      "automorphisms.json";
      probes = (isempty, x -> !isempty(x)),
    ) do r, expected
      actual = fusion_ring_automorphisms(r)

      check_equal(
        length(actual),
        length(expected) + 1,
        "automorphism count disagrees; oracle excludes identity",
      )

      check_equal(
        length(unique(actual)),
        length(actual),
        "duplicate automorphisms",
      )

      check_true(
        collect(1:rank(r)) in actual,
        "identity automorphism is missing",
      )

      N = multiplication_table(r)

      for p in actual
        valid = _prop_valid_permutation(p, rank(r))
        check_true(valid, "invalid automorphism witness")

        if valid
          check_equal_tensor(
            N[p, p, p],
            N,
            "false automorphism witness",
          )
        end
      end
    end

    _prop_reference(
      "adjoint_fusion_ring.json",
    ) do r, expected
      S, child = adjoint_fusion_ring(
        r;
        force_compute = true,
        represent_by_known = false,
      )

      _prop_check_subring(r, S, child, "reference adjoint")

      check_equal(
        length(S),
        length(expected[1]),
        "adjoint rank disagrees",
      )

      check_true(
        _prop_equivalent(
          child,
          ring_from_anyonica_code(expected[2]; fresh = true),
        ),
        "adjoint ring type disagrees",
      )
    end

    _prop_reference(
      "upper_central_series.json";
      probes = (x -> length(x) >= 3,),
    ) do r, expected
      actual = upper_central_series(
        r;
        force_compute = true,
        represent_by_known = false,
      )

      check_equal(
        length(actual),
        length(expected),
        "central-series length disagrees",
      )

      parent = r

      for (S, term) in actual
        _prop_check_subring(
          parent,
          S,
          term,
          "reference central-series embedding",
        )
        parent = term
      end

      for i in 1:min(length(actual), length(expected))
        check_true(
          _prop_equivalent(
            actual[i][2],
            ring_from_anyonica_code(
              expected[i][2];
              fresh = true,
            ),
          ),
          "central-series ring $i disagrees",
        )
      end
    end

    _prop_reference(
      "adjoint_irreps.json";
      probes = (x -> length(x) > 1,),
    ) do r, expected
      actual = adjoint_irreps(
        r;
        force_compute = true,
        represent_by_known = false,
      )

      check_equal(
        sort(length.(actual)),
        sort(length.(expected)),
        "adjoint-component sizes disagree",
      )

      check_equal(
        sort(vcat(actual...)),
        collect(1:rank(r)),
        "components do not partition the basis",
      )
    end

    _prop_reference(
      "universal_grading.json";
      probes = (
        x -> length(x[2]) == 4 && json_int(x[2][1]) > 1,
      ),
    ) do r, expected
      grades, quotient = universal_grading(
        r;
        force_compute = true,
        represent_by_known = false,
      )

      _prop_check_grading(
        r,
        grades,
        quotient,
        "reference universal grading",
      )

      reference_grades = json_int_vector(expected[1])

      actual_sizes = sort([
        count(==(g), grades)
        for g in unique(grades)
      ])

      reference_sizes = sort([
        count(==(g), reference_grades)
        for g in unique(reference_grades)
      ])

      check_equal(
        actual_sizes,
        reference_sizes,
        "grading fiber sizes disagree",
      )

      check_true(
        _prop_equivalent(
          quotient,
          ring_from_anyonica_code(expected[2]; fresh = true),
        ),
        "grading quotient disagrees",
      )
    end

    _prop_reference(
      "tensor_product_decompositions.json";
      probes = (isempty, x -> !isempty(x)),
    ) do r, expected
      actual = decompositions(
        r;
        force_compute = true,
        represent_by_known = false,
      )

      expected_rings = [
        [
          ring_from_anyonica_code(code; fresh = true)
          for code in dec
        ]
        for dec in expected
      ]

      for dec in actual
        check_true(
          all(factor -> rank(factor) >= 2, dec),
          "trivial factor in nontrivial decomposition",
        )

        check_equal(
          prod(rank.(dec)),
          rank(r),
          "factor ranks do not multiply to the parent rank",
        )
      end

      check_true(
        _prop_same_multiset(
          actual,
          expected_rings;
          equal = (a, b) -> _prop_same_multiset(a, b),
        ),
        "tensor-factor types or multiplicities disagree",
      )
    end
  end
end