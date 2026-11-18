using Test
using LinearAlgebra: kron
using Oscar: QQBar, order

# Independent tensor expectation using Kronecker products of fusion matrices.
function _ops_tensor_table(A, B)
  r, s = size(A, 1), size(B, 1)
  expected = zeros(Int, r*s, r*s, r*s)
  for a in 1:r, b in 1:s
    expected[(a-1)*s+b, :, :] = kron(A[a, :, :], B[b, :, :])
  end
  return expected
end

function _ops_check_equivalent(a, b, message)
  candidates = which_permutation(a, b; force_compute = true)
  valid = candidates !== nothing && any(
    p -> length(p) == rank(a) && sort(p) == collect(1:rank(a)) && p[1] == 1 &&
      multiplication_table(a)[p,p,p] == multiplication_table(b), candidates)
  check_true(valid, message)
  return nothing
end

@testset "Operations" begin
  maybe_testset("fusion products retain multiplicities"; level = "basic") do
    for f in small_ring_fixtures()
      @testset "$(f.name)" begin
        r, N = fixture_ring(f), f.expected.mt
        for a in 1:rank(r), b in 1:rank(r)
          expected = Dict(c => N[a,b,c] for c in 1:rank(r) if N[a,b,c] != 0)
          check_equal(FusionRings.fusion_product(r,a,b), expected, "wrong product ($a,$b)")
          check_equal(FusionRings.fusion_outcomes(r,a,b), sort(collect(keys(expected))), "wrong product support")
        end
      end
    end
    r = fixture_ring(fixture_rank2_m2())
    check_equal(FusionRings.fusion_product(r,2,2), Dict(1=>1,2=>2), "coefficient two was lost")
  end

  maybe_testset("permutation changes table and labels"; level = "basic") do
    for f in (fixture_unit(), fixture_z3(), fixture_ising(), fixture_s3())
      @testset "$(f.name)" begin
        r = fixture_ring(f)
        p = vcat(1, reverse(collect(2:rank(r))))
        expected = f.expected.mt[p,p,p]
        check_equal_tensor(FusionRings.permute_mult_tab(f.input.mt,p), expected, "wrong tensor permutation")
        result = permute(p,r)
        check_equal_tensor(multiplication_table(result), expected, "wrong ring permutation")
        check_equal(labels(result), f.expected.labels[p], "labels did not follow the basis")
        restored = permute(invperm(p),result)
        check_equal_tensor(multiplication_table(restored), f.expected.mt, "inverse permutation did not restore table")
        check_equal(labels(restored), f.expected.labels, "inverse permutation did not restore labels")
        check_equal_tensor(multiplication_table(r), f.expected.mt, "permutation mutated its input table")
        check_equal(labels(r), f.expected.labels, "permutation mutated its input labels")
      end
    end
  end

  maybe_testset("empty and single-factor tensor products"; level = "basic") do
    expected = fixture_unit().expected.mt
    check_equal_tensor(multiplication_table(tensor_product()), expected, "empty varargs product is not the unit")
    check_equal_tensor(multiplication_table(tensor_product(FusionRing[])), expected, "empty vector product is not the unit")
    r = fixture_ring(fixture_fibonacci())
    check_equal_tensor(multiplication_table(tensor_product(r)), fixture_fibonacci().expected.mt, "single-factor varargs changed the ring")
    check_equal_tensor(multiplication_table(tensor_product(FusionRing[r])), fixture_fibonacci().expected.mt, "single-factor vector changed the ring")
  end

  maybe_testset("invalid permutations"; level = "intermediate") do
    r = fixture_ring(fixture_ising())
    for p in ([1,1,3], [1,2], [2,1,3], [1,2,4])
      @testset "$p" begin
        check_throws(() -> permute(p,r), ArgumentError, "invalid ring permutation accepted")
        check_throws(() -> FusionRings.permute_mult_tab(multiplication_table(r),p), ArgumentError,
          "invalid tensor permutation accepted")
      end
    end
    check_throws(() -> FusionRings.permute_mult_tab(zeros(Int,2,2,3),[1,2]), ArgumentError,
      "noncubic tensor accepted"; message_contains = "cubic")
  end

  maybe_testset("permutations reindex cached data"; level = "intermediate") do
    f = fixture_ising()
    known_z2 = fawc(2,1,0,1)
    r = fusion_ring(copy(f.input.mt); labels = copy(f.input.labels),
      characters = copy(f.expected.characters),
      frobenius_perron_dimensions = copy(f.expected.fpdims),
      formal_codegrees = copy(f.expected.formal_codegrees),
      sub_fusion_rings = [([1,3],uuid(known_z2))],
      all_gradings = [([1,2,1],uuid(known_z2))])
    p = [1,3,2]
    result = permute(p,r)
    check_equal(result.characters, r.characters[:,p], "character columns did not follow the basis")
    check_equal(result.frobenius_perron_dimensions, r.frobenius_perron_dimensions[p], "dimension cache was not reordered")
    check_equal(result.formal_codegrees, r.formal_codegrees, "basis-independent codegrees changed")
    check_equal(sort(result.sub_fusion_rings[1][1]), [1,2], "subring embedding was not remapped")
    check_equal(result.sub_fusion_rings[1][2], uuid(known_z2), "subring UUID changed")
    check_equal(result.all_gradings[1][1], [1,1,2], "grading values were not reordered")
    check_equal(r.sub_fusion_rings[1][1], [1,3], "input subring cache was mutated")
    # A non-involutive permutation distinguishes a permutation from its inverse.
    f4 = fixture_v4()
    p4 = [1,3,4,2]
    r4 = fixture_ring(f4)
    check_equal(labels(permute(p4,r4)), f4.expected.labels[p4], "permutation direction is reversed")

    # A subring is a list of old indices: remap it with invperm(p),
    # whereas a grading is an indexed value vector: reorder it with p.
    z4 = fixture_z4()
    cached = fusion_ring(copy(z4.input.mt);
      labels = copy(z4.input.labels),
      sub_fusion_rings = [([1,3],uuid(known_z2))],
      all_gradings = [([1,2,1,2],uuid(known_z2))])
    moved = permute(p4,cached)
    check_equal(sort(moved.sub_fusion_rings[1][1]), [1,2], "three-cycle used the wrong embedding direction")
    check_equal(moved.all_gradings[1][1], [1,1,2,2], "three-cycle used the wrong grading direction")
    check_equal(cached.all_gradings[1][1], [1,2,1,2], "input grading cache was mutated")
  end

  maybe_testset("restriction checks parent closure"; level = "intermediate") do
    f = fixture_z4()
    r = fixture_ring(f)
    sub = restrict_subring(r,[3,1])
    check_equal_tensor(multiplication_table(sub), fixture_z2().expected.mt, "wrong Z4 subring")
    check_equal(labels(sub), f.expected.labels[[1,3]], "restriction lost selected labels")
    check_equal_tensor(multiplication_table(restrict_subring(r,[1,2,3,4])), f.expected.mt, "full restriction changed the ring")
    ising = fixture_ring(fixture_ising())
    check_throws(() -> restrict_subring(ising,[1,2]), ErrorException,
      "restriction silently discarded a fusion outcome"; message_contains = "not fusion-closed")
    check_throws(() -> restrict_subring(r,[1,1]), ErrorException,
      "restriction accepted duplicate indices"; message_contains = "only once")
  end

  maybe_testset("sorting and forced dimension computation"; level = "intermediate") do
    f = fixture_ising()
    r = fixture_ring(f)
    for by in ("fpdims","sd_conj"), direction in (:increasing,:decreasing)
      @testset "$by $direction" begin
        expected_p = direction == :increasing ? [1,3,2] : [1,2,3]
        result = sort(r; by = by, order = direction, force_compute = true)
        check_equal_tensor(multiplication_table(result), f.expected.mt[expected_p,expected_p,expected_p], "sort used the wrong basis order")
        check_equal(labels(result), f.expected.labels[expected_p], "sort did not reorder labels")
      end
    end
    z4 = fixture_ring(fixture_z4())
    check_equal(FusionRings.perm_vec_sd_conj(z4), [1,3,2,4], "self-duals and dual pairs were not grouped")
    stale = fusion_ring(copy(f.input.mt); frobenius_perron_dimensions = QQBar.([1,1,1]))
    check_equal(FusionRings.perm_vec_qd(stale), [1,2,3], "default sorting ignored its cache")
    check_equal(FusionRings.perm_vec_qd(stale; force_compute = true), [1,3,2], "forced sorting reused stale dimensions")
    check_equal(FusionRings.perm_vec_sd_conj(stale; force_compute = true), [1,3,2], "forced dual sorting reused stale dimensions")
    check_throws(() -> sort(r; by = "unknown"), ArgumentError, "unknown sorting mode accepted")
    check_throws(() -> sort(r; order = :sideways), ArgumentError, "unknown sorting order accepted")
  end

  maybe_testset("tensor tables, labels, and dimensions"; level = "intermediate") do
    for (a,b) in ((fixture_z2(),fixture_z3()), (fixture_fibonacci(),fixture_ising()),
                   (fixture_s3(),fixture_z2()), (fixture_rank2_m2(),fixture_z2()))
      @testset "$(a.name) tensor $(b.name)" begin
        left, right = fixture_ring(a), fixture_ring(b)
        result = tensor_product(left,right)
        check_equal_tensor(multiplication_table(result), _ops_tensor_table(a.expected.mt,b.expected.mt), "tensor product has wrong coefficients")
        check_equal(labels(result), [x*"⊗"*y for x in a.expected.labels for y in b.expected.labels], "wrong tensor basis labels")
        ds = [x*y for x in a.expected.fpdims for y in b.expected.fpdims]
        check_equal(fpdims(result), ds, "wrong stored tensor dimensions")
        check_equal(fpdims(result; force_compute = true), ds, "wrong computed tensor dimensions")
        check_equal(fpdim(result; force_compute = true), a.expected.fpdim*b.expected.fpdim, "total tensor dimension did not multiply")
        check_equal_tensor(multiplication_table(left), a.expected.mt, "left input changed")
        check_equal_tensor(multiplication_table(right), b.expected.mt, "right input changed")
      end
    end
    a,b,c = fixture_ring(fixture_z2()),fixture_ring(fixture_fibonacci()),fixture_ring(fixture_z3())
    expected = _ops_tensor_table(_ops_tensor_table(fixture_z2().expected.mt,fixture_fibonacci().expected.mt),fixture_z3().expected.mt)
    for result in (tensor_product(a,b,c), tensor_product(FusionRing[a,b,c]),
                   tensor_product(tensor_product(a,b),c), tensor_product(a,tensor_product(b,c)))
      check_equal_tensor(multiplication_table(result), expected, "tensor association or vector/varargs dispatch changed the basis")
    end
    f = fixture_ising()
    stale = fusion_ring(copy(f.input.mt); frobenius_perron_dimensions = QQBar.([1,1,1]))
    result = tensor_product(stale,fixture_ring(fixture_z2()))
    check_equal(fpdims(result), [x*y for x in f.expected.fpdims for y in fixture_z2().expected.fpdims], "tensor product propagated stale input dimensions")
  end

  maybe_testset("group conversion preserves group structure"; level = "intermediate") do
    for (f, expected_orders) in (
      (fixture_unit(), [1]),
      (fixture_z4(), [1,2,4,4]),
      (fixture_v4(), [1,2,2,2]),
      (fixture_s3(), [1,2,2,2,3,3]),
    )
      @testset "$(f.name)" begin
        group = to_group(fixture_ring(f))
        check_equal(Int(order(group)), f.expected.rank, "converted group has the wrong order")
        check_equal(sort([Int(order(g)) for g in group]), expected_orders,
          "converted group has the wrong element orders")
      end
    end
    check_throws(() -> to_group(fixture_ring(fixture_fibonacci())), ArgumentError,
      "non-group ring accepted by group conversion"; message_contains = "group fusion ring")
  end

  maybe_testset("catalogue replacement preserves basis order"; level = "intermediate") do
    f = fixture_ising()
    p = [1,3,2]
    expected = f.expected.mt[p,p,p]
    r = fusion_ring(copy(expected); frobenius_perron_dimensions = QQBar.([1,1,1]))
    kept = replace_by_known(r; keep_order = true, safe_return = false)
    check_true(kept isa FusionRing, "known Ising ring was not recognized")
    if kept isa FusionRing
      check_equal_tensor(multiplication_table(kept), expected, "keep_order changed the input basis")
    end
    canonical = replace_by_known(r; keep_order = false, safe_return = false)
    check_true(canonical isa FusionRing, "canonical replacement was not found")
    canonical isa FusionRing && _ops_check_equivalent(canonical,r,"replacement is not equivalent to input")
    # Choose an absent rank-two catalogue bucket without assuming a maximum multiplicity.
    existing = [key[2] for key in keys(FusionRings.frd) if key[1] == 2]
    m = isempty(existing) ? 1 : maximum(existing)+1
    mt = copy(fixture_rank2_m2().input.mt)
    mt[2,2,2] = m
    unknown = fusion_ring(mt)
    check_true(replace_by_known(unknown; safe_return = true) === unknown, "safe miss did not return the original ring")
    check_equal(replace_by_known(unknown; safe_return = false), missing, "unsafe miss did not return missing")
  end

  maybe_testset("Anyonica operation tables"; level = "reference") do
    @testset "permuted_tabs.json" begin
      data = load_anyonica_data("permuted_tabs.json")
      for index in oracle_case_indices("permuted_tabs.json",data)
        @testset "case $index" begin
          code, permutations = data["Input"][index]
          r = ring_from_anyonica_code(code; fresh = true)
          check_equal(length(permutations),length(data["Output"][index]),"permutation fixture has mismatched lengths")
          for j in eachindex(permutations)
            @testset "permutation $j" begin
              p = json_int_vector(permutations[j])
              expected = json_int_3tensor(data["Output"][index][j])
              check_equal_tensor(FusionRings.permute_mult_tab(multiplication_table(r),p),expected,"tensor permutation disagrees with oracle")
              result = permute(p,r)
              check_equal_tensor(multiplication_table(result),expected,"ring permutation disagrees with oracle")
              check_equal(labels(result),labels(r)[p],"reference permutation lost labels")
            end
          end
        end
      end
    end
    @testset "tensor_product_tables.json" begin
      data = load_anyonica_data("tensor_product_tables.json")
      for index in oracle_case_indices("tensor_product_tables.json",data)
        @testset "case $index" begin
          codes = data["Input"][index]
          a = ring_from_anyonica_code(codes[1]; fresh = true)
          b = ring_from_anyonica_code(codes[2]; fresh = true)
          expected = json_int_3tensor(data["Output"][index])
          check_equal_tensor(_ops_tensor_table(multiplication_table(a),multiplication_table(b)),expected,"input basis or oracle tensor convention disagrees")
          check_equal_tensor(multiplication_table(tensor_product(a,b)),expected,"tensor product disagrees with oracle")
        end
      end
    end
  end
end
