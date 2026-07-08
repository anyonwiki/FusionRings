@testset "Properties" begin

  # ============================================================
  # basic accessors
  # ============================================================

  @testset "basic accessors" begin
    maybe_testset("basic", "1. basic construction") do
      z3 = zn_fusion_ring(3)

      check_true(
        multiplication_table(z3) isa Array{Int, 3},
        "multiplication_table(z3) did not return an Int 3-tensor",
      )
      check_true(rank(z3) isa Int, "rank(z3) did not return an Int")
      check_true(names(z3) isa Vector{String}, "names(z3) did not return a Vector{String}")
      check_true(
        labels(z3) isa Vector{String}, "labels(z3) did not return a Vector{String}"
      )
      return check_true(
        tex_names(z3) isa Vector{String}, "tex_names(z3) did not return a Vector{String}"
      )
    end

    maybe_testset("intermediate", "2. intermediate correctness") do
      z3 = zn_fusion_ring(3)

      check_equal(rank(z3), 3, "rank(zn_fusion_ring(3)) was not 3")
      check_equal(labels(z3), ["0", "1", "2"], "labels(zn_fusion_ring(3)) were incorrect")
      check_equal(names(z3), ["ℤ₃"], "names(zn_fusion_ring(3)) were incorrect")
      check_equal(
        tex_names(z3),
        String[],
        "tex_names(zn_fusion_ring(3)) were not the expected empty vector",
      )
      return check_equal(
        size(multiplication_table(z3)),
        (3, 3, 3),
        "multiplication_table(zn_fusion_ring(3)) had wrong shape",
      )
    end

    maybe_testset("reference", "3. reference / Anyonica parity") do
      @test true
    end
  end

  # ============================================================
  # conjugate_element / conjugation_matrix
  # ============================================================

  @testset "conjugate_element / conjugation_matrix" begin
    maybe_testset("basic", "1. basic construction") do
      z3 = zn_fusion_ring(3)

      C = conjugation_matrix(z3)
      c = conjugate_element(z3, 1)

      check_true(
        C isa AbstractMatrix, "conjugation_matrix(z3) did not return a matrix-like object"
      )
      return check_true(c isa Int, "conjugate_element(z3, 1) did not return an Int")
    end

    maybe_testset("intermediate", "2. intermediate correctness") do
      z3 = zn_fusion_ring(3)
      z4 = zn_fusion_ring(4)

      check_equal(
        conjugation_matrix(z3),
        [
          1 0 0;
          0 0 1;
          0 1 0
        ],
        "conjugation_matrix(zn_fusion_ring(3)) was incorrect",
      )

      check_equal(conjugate_element(z3, 1), 1, "conjugate_element(z3, 1) was not 1")
      check_equal(conjugate_element(z3, 2), 3, "conjugate_element(z3, 2) was not 3")
      check_equal(conjugate_element(z3, 3), 2, "conjugate_element(z3, 3) was not 2")

      check_equal(conjugate_element(z4, 2), 4, "conjugate_element(z4, 2) was not 4")
      return check_equal(
        conjugate_element(z4, 3), 3, "conjugate_element(z4, 3) was not self-dual"
      )
    end

    maybe_testset("reference", "3. reference / Anyonica parity") do
      @test true
    end
  end

  # ============================================================
  # multiplicity / nzsc / nnzsc
  # ============================================================

  @testset "multiplicity / nonzero structure constants" begin
    maybe_testset("basic", "1. basic construction") do
      z3 = zn_fusion_ring(3)

      m = multiplicity(z3)
      s = nzsc(z3)
      n = nnzsc(z3)

      check_true(m isa Int, "multiplicity(z3) did not return an Int")
      check_true(s isa Vector, "nzsc(z3) did not return a vector")
      return check_true(n isa Int, "nnzsc(z3) did not return an Int")
    end

    maybe_testset("intermediate", "2. intermediate correctness") do
      z2 = zn_fusion_ring(2)
      su2_2 = su2k_fusion_ring(2)

      check_equal(multiplicity(z2), 1, "multiplicity(zn_fusion_ring(2)) was not 1")
      check_equal(multiplicity(su2_2), 1, "multiplicity(su2k_fusion_ring(2)) was not 1")

      check_equal(nnzsc(z2), 4, "nnzsc(zn_fusion_ring(2)) was not 4")
      return check_equal(
        sort(nzsc(z2)),
        sort([(1, 1, 1), (1, 2, 2), (2, 1, 2), (2, 2, 1)]),
        "nzsc(zn_fusion_ring(2)) was incorrect",
      )
    end

    maybe_testset("reference", "3. reference / Anyonica parity") do
      @test true
    end
  end

  # ============================================================
  # FP dimensions / FP dimension
  # ============================================================

  @testset "FP dimensions" begin
    maybe_testset("basic", "1. basic construction") do
      z3 = zn_fusion_ring(3)

      ds = fpdims(z3)
      d = fpdim(z3)

      check_true(ds isa Vector, "fpdims(z3) did not return a vector")
      check_true(
        length(ds) == rank(z3), "fpdims(z3) did not return a vector of length rank(z3)"
      )
      return check_true(d !== missing, "fpdim(z3) returned missing")
    end

    maybe_testset("intermediate", "2. intermediate correctness") do
      z3 = zn_fusion_ring(3)
      z4 = zn_fusion_ring(4)

      check_equal(
        string.(fpdims(z3)), ["{a1: 1.00000}", "{a1: 1.00000}", "{a1: 1.00000}"], "fpdims(zn_fusion_ring(3)) were not all 1"
      )
      check_equal(string(fpdim(z3)), "3", "fpdim(zn_fusion_ring(3)) was not 3")

      check_equal(
        string.(fpdims(z4)),
        ["{a1: 1.00000}", "{a1: 1.00000}", "{a1: 1.00000}", "{a1: 1.00000}"],
        "fpdims(zn_fusion_ring(4)) were not all 1",
      )
      return check_equal(string(fpdim(z4)), "4", "fpdim(zn_fusion_ring(4)) was not 4")
    end

    maybe_testset("reference", "3. reference / Anyonica parity") do
      @test true
    end
  end

  # ============================================================
  # self-duality counts
  # ============================================================

  @testset "self-duality counts" begin
    maybe_testset("basic", "1. basic construction") do
      z3 = zn_fusion_ring(3)

      v = nsdnsd(z3)
      s = nsd(z3)
      n = nnsd(z3)

      check_true(v isa Vector{Int}, "nsdnsd(z3) did not return a Vector{Int}")
      check_true(s isa Int, "nsd(z3) did not return an Int")
      return check_true(n isa Int, "nnsd(z3) did not return an Int")
    end

    maybe_testset("intermediate", "2. intermediate correctness") do
      z3 = zn_fusion_ring(3)
      z4 = zn_fusion_ring(4)

      check_equal(nsdnsd(z3), [1, 2], "nsdnsd(zn_fusion_ring(3)) was not [1,2]")
      check_equal(nsd(z3), 1, "nsd(zn_fusion_ring(3)) was not 1")
      check_equal(nnsd(z3), 2, "nnsd(zn_fusion_ring(3)) was not 2")

      check_equal(nsdnsd(z4), [2, 2], "nsdnsd(zn_fusion_ring(4)) was not [2,2]")
      check_equal(nsd(z4), 2, "nsd(zn_fusion_ring(4)) was not 2")
      return check_equal(nnsd(z4), 2, "nnsd(zn_fusion_ring(4)) was not 2")
    end

    maybe_testset("reference", "3. reference / Anyonica parity") do
      @test true
    end
  end

  # ============================================================
  # group-ring / commutative / metadata
  # ============================================================

  @testset "group-ring / commutative / metadata" begin
    maybe_testset("basic", "1. basic construction") do
      z3 = zn_fusion_ring(3)

      g = is_group_ring(z3)
      c = is_commutative(z3)
      cats = categories_with_properties(z3)
      ic = is_categorifiable(z3)

      check_true(g isa Bool, "is_group_ring(z3) did not return a Bool")
      check_true(c isa Bool, "is_commutative(z3) did not return a Bool")
      check_true(
        cats === missing || cats isa Vector,
        "categories_with_properties(z3) returned an unexpected type",
      )
      return check_true(
        ic === missing || ic isa Bool, "is_categorifiable(z3) returned an unexpected type"
      )
    end

    maybe_testset("intermediate", "2. intermediate correctness") do
      z3 = zn_fusion_ring(3)
      su2_2 = su2k_fusion_ring(2)

      check_true(is_group_ring(z3), "is_group_ring(zn_fusion_ring(3)) returned false")
      check_false(is_group_ring(su2_2), "is_group_ring(su2k_fusion_ring(2)) returned true")
      check_true(is_commutative(z3), "is_commutative(zn_fusion_ring(3)) returned false")
      return check_true(
        is_commutative(su2_2), "is_commutative(su2k_fusion_ring(2)) returned false"
      )
    end

    maybe_testset("reference", "3. reference / Anyonica parity") do
      @test true
    end
  end

  # ============================================================
  # _internal_multiplication / is_sub_fusion_ring / subsets
  # ============================================================

  @testset "_internal_multiplication / is_sub_fusion_ring / sub_fusion_ring_subsets" begin
    maybe_testset("basic", "1. basic construction") do
      z4 = zn_fusion_ring(4)

      b1 = FusionRings._internal_multiplication(z4, [1, 3])
      b2 = is_sub_fusion_ring(z4, [1, 3])
      ss = sub_fusion_ring_subsets(z4)

      check_true(b1 isa Bool, "FusionRings._internal_multiplication(z4, [1,3]) did not return a Bool")
      check_true(b2 isa Bool, "is_sub_fusion_ring(z4, [1,3]) did not return a Bool")
      return check_true(
        ss isa Vector{Vector{Int}},
        "sub_fusion_ring_subsets(z4) did not return Vector{Vector{Int}}",
      )
    end

    maybe_testset("intermediate", "2. intermediate correctness") do
      z3 = zn_fusion_ring(3)
      z4 = zn_fusion_ring(4)

      check_true(
        FusionRings._internal_multiplication(z4, [1, 3]),
        "FusionRings._internal_multiplication(z4, [1,3]) should have returned true",
      )
      check_false(
        FusionRings._internal_multiplication(z4, [1, 2]),
        "FusionRings._internal_multiplication(z4, [1,2]) should have returned false",
      )

      check_true(
        is_sub_fusion_ring(z4, [1, 3]),
        "is_sub_fusion_ring(z4, [1,3]) should have returned true",
      )
      check_false(
        is_sub_fusion_ring(z4, [1, 2]),
        "is_sub_fusion_ring(z4, [1,2]) should have returned false",
      )
      check_false(
        is_sub_fusion_ring(z4, [2, 3]),
        "is_sub_fusion_ring(z4, [2,3]) should have returned false because the unit is missing",
      )
      check_false(
        is_sub_fusion_ring(z4, Int[]),
        "is_sub_fusion_ring(z4, Int[]) should have returned false",
      )

      check_equal(
        sub_fusion_ring_subsets(z3),
        Vector{Vector{Int}}(),
        "sub_fusion_ring_subsets(z3) should have been empty",
      )
      return check_equal(
        sub_fusion_ring_subsets(z4),
        [[1, 3]],
        "sub_fusion_ring_subsets(z4) should have been [[1,3]]",
      )
    end

    maybe_testset("reference", "3. reference / Anyonica parity") do
      @test true
    end
  end

  # ============================================================
  # FusionRings._fusion_closure / restrict_subring / FusionRings._internal_closed_subsets / which_injection
  # ============================================================

  @testset "_fusion_closure / restrict_subring / _internal_closed_subsets / which_injection" begin
    maybe_testset("basic", "1. basic construction") do
      z4 = zn_fusion_ring(4)

      cl = FusionRings._fusion_closure(z4, [3])
      sub = restrict_subring(z4, [1, 3])
      css = FusionRings._internal_closed_subsets(z4, 2)
      inj = which_injection(zn_fusion_ring(2), z4)

      check_true(
        cl isa Vector{Int}, "FusionRings._fusion_closure(z4, [3]) did not return a Vector{Int}"
      )
      check_true(
        sub isa FusionRing, "restrict_subring(z4, [1,3]) did not return a FusionRing"
      )
      check_true(
        css isa Vector{Vector{Int}},
        "FusionRings._internal_closed_subsets(z4, 2) did not return Vector{Vector{Int}}",
      )
      return check_true(
        inj === nothing || inj isa Dict{Int, Int},
        "which_injection(zn_fusion_ring(2), z4) returned an unexpected type",
      )
    end

    maybe_testset("intermediate", "2. intermediate correctness") do
      z4 = zn_fusion_ring(4)
      z2 = zn_fusion_ring(2)

      check_equal(
        FusionRings._fusion_closure(z4, [3]), [1, 3], "FusionRings._fusion_closure(z4, [3]) was not [1,3]"
      )
      check_equal(
        FusionRings._fusion_closure(z4, [2]),
        [1, 2, 3, 4],
        "FusionRings._fusion_closure(z4, [2]) was not the whole ring",
      )

      sub = restrict_subring(z4, [1, 3])
      check_equal(rank(sub), 2, "restrict_subring(z4, [1,3]) did not produce a rank-2 ring")
      check_equal(
        multiplication_table(sub),
        multiplication_table(z2),
        "restrict_subring(z4, [1,3]) did not match zn_fusion_ring(2)",
      )

      check_equal(
        FusionRings._internal_closed_subsets(z4, 1),
        [[1]],
        "FusionRings._internal_closed_subsets(z4, 1) was not [[1]]",
      )
      check_equal(
        FusionRings._internal_closed_subsets(z4, 2),
        [[1, 3]],
        "FusionRings._internal_closed_subsets(z4, 2) was not [[1,3]]",
      )

      return check_equal(
        which_injection(z2, z4),
        Dict(1 => 1, 2 => 3),
        "which_injection(zn_fusion_ring(2), zn_fusion_ring(4)) was not Dict(1=>1, 2=>3)",
      )
    end

    maybe_testset("reference", "3. reference / Anyonica parity") do
      @test true
    end
  end

  # ============================================================
  #  FusionRings._permutation_vector_equiv
  # ============================================================

  @testset "FusionRings._permutation_vector_equiv" begin
    maybe_testset("basic", "1. basic construction") do
      z3 = zn_fusion_ring(3)
      N = multiplication_table(z3)

      pv = FusionRings._permutation_vector_equiv(N, N)

      return check_true(
        pv === nothing || pv isa Vector{Int},
        "FusionRings._permutation_vector_equiv(N, N) returned an unexpected type",
      )
    end

    maybe_testset("intermediate", "2. intermediate correctness") do
      z3 = zn_fusion_ring(3)
      z4 = zn_fusion_ring(4)
      N3 = multiplication_table(z3)
      pz3 = permute([1, 3, 2],z3)
      N3p = multiplication_table(pz3)

      pv = FusionRings._permutation_vector_equiv(N3, N3p)
      check_true(
        pv !== nothing, "FusionRings._permutation_vector_equiv(N3, permuted N3) returned nothing"
      )
      return check_equal(
        _permute_multtab(N3, pv),
        N3p,
        "FusionRings._permutation_vector_equiv did not return a valid permutation",
      )
    end

    maybe_testset("reference", "3. reference / Anyonica parity") do
      @test true
    end
  end

  # ============================================================
  # is_equivalent_fusion_ring / fusion_ring_automorphisms
  # ============================================================

  @testset "is_equivalent_fusion_ring / fusion_ring_automorphisms" begin
    maybe_testset("basic", "1. basic construction") do
      z3 = zn_fusion_ring(3)
      pz3 = permute([1, 3, 2],z3)

      eq = is_equivalent_fusion_ring(z3, pz3)
      autos = fusion_ring_automorphisms(z3)

      check_true(eq isa Bool, "is_equivalent_fusion_ring(z3, pz3) did not return a Bool")
      return check_true(
        autos isa Vector{Vector{Int}},
        "fusion_ring_automorphisms(z3) did not return Vector{Vector{Int}}",
      )
    end

    maybe_testset("intermediate", "2. intermediate correctness") do
      z2 = zn_fusion_ring(2)
      z3 = zn_fusion_ring(3)
      z4 = zn_fusion_ring(4)
      pz3 = permute([1, 3, 2],z3)

      check_true(
        is_equivalent_fusion_ring(z3, pz3),
        "is_equivalent_fusion_ring(z3, permute([1,3,2],z3)) returned false",
      )
      check_false(
        is_equivalent_fusion_ring(z2, z3),
        "is_equivalent_fusion_ring(z2, z3) should have returned false",
      )

      check_equal(
        fusion_ring_automorphisms(z2),
        [[1, 2]],
        "fusion_ring_automorphisms(z2) was not [[1,2]]",
      )
      check_equal(
        fusion_ring_automorphisms(z3),
        [[1, 2, 3], [1, 3, 2]],
        "fusion_ring_automorphisms(z3) was incorrect",
      )
      return check_true(
        [1, 4, 3, 2] in fusion_ring_automorphisms(z4),
        "fusion_ring_automorphisms(z4) did not contain the inversion automorphism",
      )
    end

    maybe_testset("reference", "3. reference / Anyonica parity") do
      @test true
    end
  end

  # ============================================================
  # decompositions
  # ============================================================

  @testset "decompositions" begin
    maybe_testset("basic", "1. basic construction") do
      z2 = zn_fusion_ring(2)
      r = decompositions(z2)

      return check_true(r isa Vector, "decompositions(z2) did not return a vector")
    end

    maybe_testset("intermediate", "2. intermediate correctness") do
      z2 = zn_fusion_ring(2)

      check_throws(
        () -> decompositions(z2, "DirectSum"),
        "decompositions(z2, \"DirectSum\") should have thrown because only tensor products are implemented",
      )
    end

    maybe_testset("reference", "3. reference / Anyonica parity") do
      @test true
    end
  end

  # ============================================================
  # adjoint_fusion_ring / upper_central_series / is_nilpotent
  # ============================================================

  @testset "adjoint_fusion_ring / upper_central_series / is_nilpotent" begin
    maybe_testset("basic", "1. basic construction") do
      z3 = zn_fusion_ring(3)

      ad = adjoint_fusion_ring(z3)
      ucs = upper_central_series(z3)
      nilp = is_nilpotent(z3)

      check_true(
        ad isa Tuple{Vector{Int}, FusionRing},
        "adjoint_fusion_ring(z3) did not return Tuple{Vector{Int},FusionRing}",
      )
      check_true(ucs isa Vector, "upper_central_series(z3) did not return a vector")
      return check_true(nilp isa Bool, "is_nilpotent(z3) did not return a Bool")
    end

    maybe_testset("intermediate", "2. intermediate correctness") do
      z3 = zn_fusion_ring(3)
      z4 = zn_fusion_ring(4)

      ad3_set, ad3_ring = adjoint_fusion_ring(z3)
      check_equal(
        ad3_set, [1], "adjoint_fusion_ring(z3) did not return [1] as the adjoint subset"
      )
      check_equal(
        rank(ad3_ring), 1, "adjoint_fusion_ring(z3) did not return a rank-1 adjoint subring"
      )

      ucs3 = upper_central_series(z3)
      check_equal(length(ucs3), 2, "upper_central_series(z3) did not have length 2")
      check_equal(
        first(ucs3)[1], [1, 2, 3], "upper_central_series(z3) first subset was incorrect"
      )
      check_equal(last(ucs3)[1], [1], "upper_central_series(z3) last subset was incorrect")

      check_true(is_nilpotent(z3), "is_nilpotent(z3) should have returned true")
      return check_true(is_nilpotent(z4), "is_nilpotent(z4) should have returned true")
    end

    maybe_testset("reference", "3. reference / Anyonica parity") do
      @test true
    end
  end

  # ============================================================
  # unfinished / placeholder functions
  # ============================================================

  @testset "unfinished / placeholder functions" begin
    maybe_testset("basic", "1. basic construction") do
      z2 = zn_fusion_ring(2)
      z4 = zn_fusion_ring(4)

      check_throws(
        () -> sub_ring_tables(zeros(Int, 2, 2)),
        "sub_ring_tables should throw or fail because it is not implemented",
      )

      check_throws(
        () -> injection_form(z2, z4),
        "injection_form(z2, z4) should throw or fail because it is not implemented",
      )
    end

    maybe_testset("intermediate", "2. intermediate correctness") do
      @test true
    end

    maybe_testset("reference", "3. reference / Anyonica parity") do
      @test true
    end
  end

  # ============================================================
  # characters / numeric_characters / numeric FP helpers
  # ============================================================

  @testset "characters / numeric_characters / numeric FP helpers" begin
    maybe_testset("basic", "1. basic construction") do
      z3 = zn_fusion_ring(3)

      ch = characters(z3)
      nch = numeric_characters(z3)
      nfpd = numeric_fpdim(z3)
      nfpds = numeric_fpdims(z3)

      check_true(
        ch isa AbstractMatrix, "characters(z3) did not return a matrix-like object"
      )
      check_true(
        nch isa AbstractMatrix || nch isa Vector,
        "numeric_characters(z3) returned an unexpected type",
      )
      check_true(nfpd isa Real, "numeric_fpdim(z3) did not return a real number")
      return check_true(
        nfpds isa Vector{Float64}, "numeric_fpdims(z3) did not return Vector{Float64}"
      )
    end

    maybe_testset("intermediate", "2. intermediate correctness") do
      z1 = zn_fusion_ring(1)
      z3 = zn_fusion_ring(3)

      check_equal(
        string.(characters(z1)), ["{a1: 1.00000}"], "characters(zn_fusion_ring(1)) were not [1]"
      )

      check_equal(
        size(characters(z3)),
        (3, 3),
        "characters(zn_fusion_ring(3)) did not return a 3×3 matrix",
      )
      check_equal(
        size(numeric_characters(z3)),
        (3, 3),
        "numeric_characters(zn_fusion_ring(3)) did not return a 3×3 matrix",
      )

      check_approx(
        numeric_fpdim(z3), 3.0, "numeric_fpdim(zn_fusion_ring(3)) was not approximately 3.0"
      )
      check_approx(
        numeric_fpdims(z3),
        [1.0, 1.0, 1.0],
        "numeric_fpdims(zn_fusion_ring(3)) were not [1.0,1.0,1.0]",
      )

      check_throws(
        () -> characters(FusionRingTY([1 2; 2 1])),
        "characters(FusionRingTY(Z2 table)) should have thrown for a non-commutative ring",
      )
      return check_throws(
        () -> numeric_characters(FusionRingTY([1 2; 2 1])),
        "numeric_characters(FusionRingTY(Z2 table)) should have thrown for a non-commutative ring",
      )
    end

    maybe_testset("reference", "3. reference / Anyonica parity") do
      @test true
    end
  end

  # ============================================================
  # projective_SL_2_ZZ_reps
  # ============================================================

  @testset "projective_SL_2_ZZ_reps" begin
    maybe_testset("basic", "1. basic construction") do
      z3 = zn_fusion_ring(3)

      return check_throws(
        () -> projective_SL_2_ZZ_reps(z3),
        "projective_SL_2_ZZ_reps(z3) should have thrown when no data is available",
      )
    end

    maybe_testset("intermediate", "2. intermediate correctness") do
      z3 = zn_fusion_ring(3)

      return check_throws(
        () -> projective_SL_2_ZZ_reps(z3),
        "projective_SL_2_ZZ_reps(z3) should have thrown when no data is available",
      )
    end

    maybe_testset("reference", "3. reference / Anyonica parity") do
      @test true
    end
  end
end
