@testset "Operations" begin

  #TODO (things i fixed:)
  #1. Uses permute(p, r), not permute(r, p).
  #2. Uses is_equivalent_fusion_ring, not is_equivalent.
  #3. Uses first_permutation_result(...) because which_permutation may return nothing or a vector/list of vectors.
  #4. Uses formal codes from the operation JSON files to reconstruct stored rings with from_anyonwiki_code.
  #5. Uses check_equal_tensor/check_mt_equal for full multiplication-table comparisons.
  #6. Qualifies internal permute_mult_tab as FusionRings.permute_mult_tab.
  # ============================================================
  # fusion_matrix / fusion_coeff / fusion_product
  # ============================================================

  @testset "fusion_matrix / fusion_coeff / fusion_product" begin
    maybe_testset("basic", "1. basic construction") do
      z2 = zn_fusion_ring(2)

      fm = fusion_matrix(z2, 1)
      fc = fusion_coeff(z2, 1, 1, 1)
      fp = fusion_product(z2, 1, 1)

      check_true(fm isa Matrix{Int}, "fusion_matrix(z2, 1) did not return a Matrix{Int}")
      check_true(fc isa Int, "fusion_coeff(z2, 1, 1, 1) did not return an Int")
      return check_true(
        fp isa Dict{Int, Int}, "fusion_product(z2, 1, 1) did not return a Dict{Int,Int}"
      )
    end

    maybe_testset("intermediate", "2. intermediate correctness") do
      z2 = zn_fusion_ring(2)
      z3 = zn_fusion_ring(3)

      check_equal(
        fusion_matrix(z2, 1), [1 0; 0 1], "fusion_matrix(z2, 1) was not the identity matrix"
      )

      check_equal(
        fusion_matrix(z2, 2),
        [0 1; 1 0],
        "fusion_matrix(z2, 2) was not the expected permutation matrix",
      )

      check_equal(fusion_coeff(z2, 1, 1, 1), 1, "fusion_coeff(z2, 1, 1, 1) was not 1")

      check_equal(fusion_coeff(z2, 2, 2, 1), 1, "fusion_coeff(z2, 2, 2, 1) was not 1")

      check_equal(fusion_coeff(z2, 2, 2, 2), 0, "fusion_coeff(z2, 2, 2, 2) was not 0")

      check_equal(
        fusion_product(z2, 1, 1),
        Dict(1 => 1),
        "fusion_product(z2, 1, 1) was not Dict(1 => 1)",
      )

      check_equal(
        fusion_product(z2, 2, 2),
        Dict(1 => 1),
        "fusion_product(z2, 2, 2) was not Dict(1 => 1)",
      )

      check_equal(
        fusion_product(z3, 2, 2),
        Dict(3 => 1),
        "fusion_product(z3, 2, 2) was not Dict(3 => 1)",
      )

      return check_equal(
        fusion_product(z3, 2, 3),
        Dict(1 => 1),
        "fusion_product(z3, 2, 3) was not Dict(1 => 1)",
      )
    end

    maybe_testset("reference", "3. reference / Anyonica parity") do
      # No direct Anyonica fixture was exported for individual
      # fusion_matrix/fusion_coeff/fusion_product calls.
      #
      # These are tested indirectly through the constructor,
      # permutation, and tensor-product multiplication-table fixtures.
      @test true
    end
  end

  # ============================================================
  # fusion_outcomes / decompose
  # ============================================================

  @testset "fusion_outcomes / decompose" begin
    maybe_testset("basic", "1. basic construction") do
      z2 = zn_fusion_ring(2)

      fo = fusion_outcomes(z2, 1, 1)
      dec = decompose(z2, 1, 1)

      check_true(
        fo isa Vector{Int}, "fusion_outcomes(z2, 1, 1) did not return a Vector{Int}"
      )
      return check_true(dec isa Vector, "decompose(z2, 1, 1) did not return a vector")
    end

    maybe_testset("intermediate", "2. intermediate correctness") do
      z2 = zn_fusion_ring(2)
      z3 = zn_fusion_ring(3)
      su2_2 = su2k_fusion_ring(2)

      check_equal(fusion_outcomes(z2, 1, 1), [1], "fusion_outcomes(z2, 1, 1) was not [1]")

      check_equal(fusion_outcomes(z2, 2, 2), [1], "fusion_outcomes(z2, 2, 2) was not [1]")

      check_equal(fusion_outcomes(z3, 2, 2), [3], "fusion_outcomes(z3, 2, 2) was not [3]")

      check_equal(
        fusion_outcomes(su2_2, 2, 2), [1, 3], "fusion_outcomes(su2_2, 2, 2) was not [1,3]"
      )

      check_equal(decompose(z2, 2, 2), [(1, 1)], "decompose(z2, 2, 2) was not [(1,1)]")

      return check_equal(
        decompose(su2_2, 2, 2),
        [(1, 1), (3, 1)],
        "decompose(su2_2, 2, 2) was not [(1,1), (3,1)]",
      )
    end

    maybe_testset("reference", "3. reference / Anyonica parity") do
      # No direct Anyonica fixture was exported for fusion_outcomes/decompose.
      # The tensor-product and permutation fixtures test the same
      # multiplication-table data these functions read from.
      @test true
    end
  end

  # ============================================================
  # permute_mult_tab
  # ============================================================

  @testset "permute_mult_tab" begin
    maybe_testset("basic", "1. basic construction") do
      z2 = zn_fusion_ring(2)
      N = multiplication_table(z2)

      M = FusionRings.permute_mult_tab(N, [1, 2])

      check_true(M isa Array{Int, 3}, "permute_mult_tab did not return an Int 3-tensor")
      return check_equal(
        size(M), size(N), "permute_mult_tab did not preserve the multiplication-table shape"
      )
    end

    maybe_testset("intermediate", "2. intermediate correctness") do
      z2 = zn_fusion_ring(2)
      z3 = zn_fusion_ring(3)

      N2 = multiplication_table(z2)
      N3 = multiplication_table(z3)

      check_equal_tensor(
        FusionRings.permute_mult_tab(N2, [1, 2]),
        N2,
        "permute_mult_tab(N2, [1,2]) did not leave the table unchanged",
      )

      check_equal_tensor(
        FusionRings.permute_mult_tab(N3, [1, 3, 2]),
        N3,
        "permute_mult_tab on Z3 with permutation [1,3,2] did not preserve the table as expected",
      )

      return check_throws(
        () -> FusionRings.permute_mult_tab(N3, [2, 1, 3]),
        "permute_mult_tab accepted a permutation that did not fix the vacuum",
      )
    end

    maybe_testset("reference", "3. reference / Anyonica parity") do
      data = load_anyonica_data("permuted_tabs.json")

      for idx in oracle_case_indices("permuted_tabs.json", data)
        code = data["Input"][idx][1]
        perms = data["Input"][idx][2]
        expected_tabs = data["Output"][idx]

        r = ring_from_anyonica_code(code)
        N = multiplication_table(r)

        for j in eachindex(perms)
          p = json_int_vector(perms[j])

          actual = FusionRings.permute_mult_tab(N, p)

          check_mt_equal(
            actual,
            expected_tabs[j],
            "permute_mult_tab did not match Anyonica permuted_tabs.json case ($idx, $j)",
          )
        end
      end
    end
  end

  # ============================================================
  # permute
  # ============================================================

  @testset "permute" begin
    maybe_testset("basic", "1. basic construction") do
      z3 = zn_fusion_ring(3)
      pz3 = permute([1, 3, 2], z3)

      check_true(pz3 isa FusionRing, "permute([1,3,2], z3) did not return a FusionRing")
      check_equal(
        size(multiplication_table(pz3)),
        (rank(pz3), rank(pz3), rank(pz3)),
        "permute([1,3,2], z3) did not produce a cubic multiplication table",
      )
      return check_equal(
        length(labels(pz3)), rank(pz3), "permute([1,3,2], z3) did not preserve label count"
      )
    end

    maybe_testset("intermediate", "2. intermediate correctness") do
      z3 = zn_fusion_ring(3)
      pz3 = permute([1, 3, 2], z3)

      check_equal(
        labels(pz3),
        ["0", "2", "1"],
        "permute([1,3,2], z3) did not permute labels correctly",
      )

      check_true(
        is_equivalent_fusion_ring(z3, pz3),
        "permute([1,3,2], z3) did not produce an equivalent fusion ring",
      )

      check_throws(
        () -> permute([1, 1, 2], z3), "permute accepted a non-permutation vector"
      )

      check_throws(
        () -> permute([2, 1, 3], z3),
        "permute accepted a permutation that did not fix the vacuum",
      )

      return check_throws(
        () -> permute([1, 2], z3), "permute accepted a permutation of the wrong length"
      )
    end

    maybe_testset("reference", "3. reference / Anyonica parity") do
      data = load_anyonica_data("permuted_tabs.json")

      for idx in oracle_case_indices("permuted_tabs.json", data)
        code = data["Input"][idx][1]
        perms = data["Input"][idx][2]
        expected_tabs = data["Output"][idx]

        r = ring_from_anyonica_code(code)

        for j in eachindex(perms)
          p = json_int_vector(perms[j])

          actual = permute(p, r)

          check_mt_equal(
            actual,
            expected_tabs[j],
            "permute did not match Anyonica permuted_tabs.json case ($idx, $j)",
          )
        end
      end
    end
  end

  # ============================================================
  # perm_vec_qd / perm_vec_sd_conj
  # ============================================================

  @testset "perm_vec_qd / perm_vec_sd_conj" begin
    maybe_testset("basic", "1. basic construction") do
      z3 = zn_fusion_ring(3)
      qd_perm = perm_vec_qd(z3)
      sd_perm = perm_vec_sd_conj(z3)

      check_true(qd_perm isa Vector{Int}, "perm_vec_qd(z3) did not return a Vector{Int}")
      return check_true(
        sd_perm isa Vector{Int}, "perm_vec_sd_conj(z3) did not return a Vector{Int}"
      )
    end

    maybe_testset("intermediate", "2. intermediate correctness") do
      z3 = zn_fusion_ring(3)
      z4 = zn_fusion_ring(4)

      check_equal(
        perm_vec_qd(z3), [1, 2, 3], "perm_vec_qd(z3) was not the identity permutation"
      )

      check_equal(
        perm_vec_qd(z3; order = :decreasing),
        [1, 2, 3],
        "perm_vec_qd(z3; order=:decreasing) was not the identity permutation",
      )

      check_equal(
        perm_vec_sd_conj(z3),
        [1, 2, 3],
        "perm_vec_sd_conj(z3) was not the identity permutation",
      )

      return check_equal(
        perm_vec_sd_conj(z4),
        [1, 3, 2, 4],
        "perm_vec_sd_conj(z4) did not place the self-dual object before the conjugate pair",
      )
    end

    maybe_testset("reference", "3. reference / Anyonica parity") do
      # No direct Anyonica fixture was exported for these sorting helpers.
      @test true
    end
  end

  # ============================================================
  # tensor_product
  # ============================================================

  @testset "tensor_product" begin
    maybe_testset("basic", "1. basic construction") do
      z2 = zn_fusion_ring(2)
      z3 = zn_fusion_ring(3)
      tp = tensor_product(z2, z3)

      check_true(tp isa FusionRing, "tensor_product(z2, z3) did not return a FusionRing")
      check_equal(
        size(multiplication_table(tp)),
        (rank(tp), rank(tp), rank(tp)),
        "tensor_product(z2, z3) did not produce a cubic multiplication table",
      )
      return check_equal(
        length(labels(tp)),
        rank(tp),
        "tensor_product(z2, z3) did not produce the correct number of labels",
      )
    end

    maybe_testset("intermediate", "2. intermediate correctness") do
      z2 = zn_fusion_ring(2)
      z3 = zn_fusion_ring(3)
      tp = tensor_product(z2, z3)

      check_equal(rank(tp), 6, "tensor_product(z2, z3) did not have rank 6")

      check_equal(
        labels(tp),
        ["0⊗0", "0⊗1", "0⊗2", "1⊗0", "1⊗1", "1⊗2"],
        "tensor_product(z2, z3) did not produce the expected labels",
      )

      check_equal(
        fusion_product(tp, 1, 1),
        Dict(1 => 1),
        "tensor_product(z2, z3): vacuum × vacuum was not vacuum",
      )

      return check_equal(
        fusion_product(tp, 4, 4),
        Dict(1 => 1),
        "tensor_product(z2, z3): (1⊗0) × (1⊗0) was not 0⊗0",
      )
    end

    maybe_testset("reference", "3. reference / Anyonica parity") do
      data = load_anyonica_data("tensor_product_tables.json")

      for idx in oracle_case_indices("tensor_product_tables.json", data)
        code1 = data["Input"][idx][1]
        code2 = data["Input"][idx][2]

        r1 = ring_from_anyonica_code(code1)
        r2 = ring_from_anyonica_code(code2)

        actual = tensor_product(r1, r2)

        check_mt_equal(
          actual,
          data["Output"][idx],
          "tensor_product did not match Anyonica tensor_product_tables.json case $idx",
        )
      end
    end
  end

  # ============================================================
  # which_permutation / is_equivalent_fusion_ring
  # ============================================================

  @testset "which_permutation / is_equivalent_fusion_ring" begin
    maybe_testset("basic", "1. basic construction") do
      z3 = zn_fusion_ring(3)
      pz3 = permute([1, 3, 2], z3)

      wp = first_permutation_result(which_permutation(z3, pz3))
      eq = is_equivalent_fusion_ring(z3, pz3)

      check_true(
        wp === nothing || wp isa Vector{Int},
        "which_permutation(z3, pz3) returned neither nothing nor a permutation vector",
      )
      return check_true(
        eq isa Bool, "is_equivalent_fusion_ring(z3, pz3) did not return a Bool"
      )
    end

    maybe_testset("intermediate", "2. intermediate correctness") do
      z2 = zn_fusion_ring(2)
      z3 = zn_fusion_ring(3)
      z4 = zn_fusion_ring(4)
      pz3 = permute([1, 3, 2], z3)

      wp = first_permutation_result(which_permutation(z3, pz3))

      check_true(
        wp !== nothing, "which_permutation(z3, permute([1,3,2], z3)) returned nothing"
      )

      check_equal_tensor(
        FusionRings.permute_mult_tab(multiplication_table(z3), wp),
        multiplication_table(pz3),
        "which_permutation did not return a valid matching permutation",
      )

      check_true(
        is_equivalent_fusion_ring(z3, pz3),
        "is_equivalent_fusion_ring(z3, permute([1,3,2], z3)) returned false",
      )

      check_equal(
        which_permutation(z2, z4),
        nothing,
        "which_permutation(z2, z4) should have returned nothing for different ranks",
      )

      return check_false(
        is_equivalent_fusion_ring(z2, z4),
        "is_equivalent_fusion_ring(z2, z4) should have returned false",
      )
    end

    maybe_testset("reference", "3. reference / Anyonica parity") do
      data = load_anyonica_data("permuted_tabs.json")

      for idx in oracle_case_indices("permuted_tabs.json", data)
        code = data["Input"][idx][1]
        expected_tabs = data["Output"][idx]

        r = ring_from_anyonica_code(code)

        for j in eachindex(expected_tabs)
          target = ring_from_anyonica_mt(expected_tabs[j])

          p = first_permutation_result(which_permutation(r, target))

          check_true(
            p !== nothing,
            "which_permutation returned nothing for Anyonica permuted_tabs.json case ($idx, $j)",
          )

          check_equal_tensor(
            FusionRings.permute_mult_tab(multiplication_table(r), p),
            multiplication_table(target),
            "which_permutation returned an invalid permutation for Anyonica permuted_tabs.json case ($idx, $j)",
          )

          check_true(
            is_equivalent_fusion_ring(r, target),
            "is_equivalent_fusion_ring returned false for Anyonica permuted_tabs.json case ($idx, $j)",
          )
        end
      end
    end
  end
end
