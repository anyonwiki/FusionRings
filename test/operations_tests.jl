@testset "Operations" begin

    # ============================================================
    # fusion_matrix / fusion_coeff / fusion_product
    # ============================================================

    @testset "fusion_matrix / fusion_coeff / fusion_product" begin
        maybe_testset("basic", "1. basic construction") do
            z2 = zn_fusion_ring(2)

            fm = fusion_matrix(z2, 1)
            fc = fusion_coeff(z2, 1, 1, 1)
            fp = fusion_product(z2, 1, 1)

            check_true(fm isa Matrix{Int},
                "fusion_matrix(z2, 1) did not return a Matrix{Int}")
            check_true(fc isa Int,
                "fusion_coeff(z2, 1, 1, 1) did not return an Int")
            check_true(fp isa Dict{Int,Int},
                "fusion_product(z2, 1, 1) did not return a Dict{Int,Int}")
        end

        maybe_testset("intermediate", "2. intermediate correctness") do
            z2 = zn_fusion_ring(2)
            z3 = zn_fusion_ring(3)

            check_equal(
                fusion_matrix(z2, 1),
                [1 0; 0 1],
                "fusion_matrix(z2, 1) was not the identity matrix"
            )

            check_equal(
                fusion_matrix(z2, 2),
                [0 1; 1 0],
                "fusion_matrix(z2, 2) was not the expected permutation matrix"
            )

            check_equal(
                fusion_coeff(z2, 1, 1, 1),
                1,
                "fusion_coeff(z2, 1, 1, 1) was not 1"
            )

            check_equal(
                fusion_coeff(z2, 2, 2, 1),
                1,
                "fusion_coeff(z2, 2, 2, 1) was not 1"
            )

            check_equal(
                fusion_coeff(z2, 2, 2, 2),
                0,
                "fusion_coeff(z2, 2, 2, 2) was not 0"
            )

            check_equal(
                fusion_product(z2, 1, 1),
                Dict(1 => 1),
                "fusion_product(z2, 1, 1) was not Dict(1 => 1)"
            )

            check_equal(
                fusion_product(z2, 2, 2),
                Dict(1 => 1),
                "fusion_product(z2, 2, 2) was not Dict(1 => 1)"
            )

            check_equal(
                fusion_product(z3, 2, 2),
                Dict(3 => 1),
                "fusion_product(z3, 2, 2) was not Dict(3 => 1)"
            )

            check_equal(
                fusion_product(z3, 2, 3),
                Dict(1 => 1),
                "fusion_product(z3, 2, 3) was not Dict(1 => 1)"
            )
        end

        maybe_testset("reference", "3. reference / Anyonica parity") do
            # TODO
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

            check_true(fo isa Vector{Int},
                "fusion_outcomes(z2, 1, 1) did not return a Vector{Int}")
            check_true(dec isa Vector,
                "decompose(z2, 1, 1) did not return a vector")
        end

        maybe_testset("intermediate", "2. intermediate correctness") do
            z2 = zn_fusion_ring(2)
            z3 = zn_fusion_ring(3)
            su2_2 = su2k_fusion_ring(2)

            check_equal(
                fusion_outcomes(z2, 1, 1),
                [1],
                "fusion_outcomes(z2, 1, 1) was not [1]"
            )

            check_equal(
                fusion_outcomes(z2, 2, 2),
                [1],
                "fusion_outcomes(z2, 2, 2) was not [1]"
            )

            check_equal(
                fusion_outcomes(z3, 2, 2),
                [3],
                "fusion_outcomes(z3, 2, 2) was not [3]"
            )

            check_equal(
                fusion_outcomes(su2_2, 2, 2),
                [1, 3],
                "fusion_outcomes(su2_2, 2, 2) was not [1,3]"
            )

            check_equal(
                decompose(z2, 2, 2),
                [(1,1)],
                "decompose(z2, 2, 2) was not [(1,1)]"
            )

            check_equal(
                decompose(su2_2, 2, 2),
                [(1,1), (3,1)],
                "decompose(su2_2, 2, 2) was not [(1,1), (3,1)]"
            )
        end

        maybe_testset("reference", "3. reference / Anyonica parity") do
            # TODO
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

            M = permute_mult_tab(N, [1,2])

            check_true(M isa Array{Int,3},
                "permute_mult_tab did not return an Int 3-tensor")
            check_equal(size(M), size(N),
                "permute_mult_tab did not preserve the multiplication-table shape")
        end

        maybe_testset("intermediate", "2. intermediate correctness") do
            z2 = zn_fusion_ring(2)
            z3 = zn_fusion_ring(3)

            N2 = multiplication_table(z2)
            N3 = multiplication_table(z3)

            check_equal(
                permute_mult_tab(N2, [1,2]),
                N2,
                "permute_mult_tab(N2, [1,2]) did not leave the table unchanged"
            )

            check_equal(
                permute_mult_tab(N3, [1,3,2]),
                N3,
                "permute_mult_tab on Z3 with permutation [1,3,2] did not preserve the table as expected"
            )

            check_throws(
                () -> permute_mult_tab(N3, [2,1,3]),
                "permute_mult_tab accepted a permutation that did not fix the vacuum"
            )
        end

        maybe_testset("reference", "3. reference / Anyonica parity") do
            # TODO
            @test true
        end
    end

    # ============================================================
    # permute
    # ============================================================

    @testset "permute" begin
        maybe_testset("basic", "1. basic construction") do
            z3 = zn_fusion_ring(3)
            pz3 = permute( [1,3,2], z3)

            check_true(pz3 isa FusionRing,
                "permute([1,3,2],z3) did not return a FusionRing")
            check_equal(
                size(multiplication_table(pz3)),
                (rank(pz3), rank(pz3), rank(pz3)),
                "permute([1,3,2],z3) did not produce a cubic multiplication table"
            )
            check_equal(length(labels(pz3)), rank(pz3),
                "permute([1,3,2],z3) did not preserve label count")
        end

        maybe_testset("intermediate", "2. intermediate correctness") do
            z3 = zn_fusion_ring(3)
            pz3 = permute([1,3,2],z3)

            check_equal(
                labels(pz3),
                ["0", "2", "1"],
                "permute([1,3,2],z3) did not permute labels correctly"
            )

            check_true(
                is_equivalent(z3, pz3),
                "permute([1,3,2],z3) did not produce an equivalent fusion ring"
            )

            check_throws(
                () -> permute( [1,1,2],z3),
                "permute accepted a non-permutation vector"
            )

            check_throws(
                () -> permute([2,1,3],z3),
                "permute accepted a permutation that did not fix the vacuum"
            )

            check_throws(
                () -> permute([1,2],z3),
                "permute accepted a permutation of the wrong length"
            )
        end

        maybe_testset("reference", "3. reference / Anyonica parity") do
            # TODO
            @test true
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

            check_true(qd_perm isa Vector{Int},
                "perm_vec_qd(z3) did not return a Vector{Int}")
            check_true(sd_perm isa Vector{Int},
                "perm_vec_sd_conj(z3) did not return a Vector{Int}")
        end

        maybe_testset("intermediate", "2. intermediate correctness") do
            z3 = zn_fusion_ring(3)
            z4 = zn_fusion_ring(4)

            check_equal(
                perm_vec_qd(z3),
                [1,2,3],
                "perm_vec_qd(z3) was not the identity permutation"
            )

            check_equal(
                perm_vec_qd(z3; order=:decreasing),
                [1,2,3],
                "perm_vec_qd(z3; order=:decreasing) was not the identity permutation"
            )

            check_equal(
                perm_vec_sd_conj(z3),
                [1,2,3],
                "perm_vec_sd_conj(z3) was not the identity permutation"
            )

            check_equal(
                perm_vec_sd_conj(z4),
                [1,3,2,4],
                "perm_vec_sd_conj(z4) did not place the self-dual object before the conjugate pair"
            )
        end

        maybe_testset("reference", "3. reference / Anyonica parity") do
            # TODO
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

            check_true(tp isa FusionRing,
                "tensor_product(z2, z3) did not return a FusionRing")
            check_equal(
                size(multiplication_table(tp)),
                (rank(tp), rank(tp), rank(tp)),
                "tensor_product(z2, z3) did not produce a cubic multiplication table"
            )
            check_equal(length(labels(tp)), rank(tp),
                "tensor_product(z2, z3) did not produce the correct number of labels")
        end

        maybe_testset("intermediate", "2. intermediate correctness") do
            z2 = zn_fusion_ring(2)
            z3 = zn_fusion_ring(3)
            tp = tensor_product(z2, z3)

            check_equal(
                rank(tp),
                6,
                "tensor_product(z2, z3) did not have rank 6"
            )

            check_equal(
                labels(tp),
                ["0⊗0", "0⊗1", "0⊗2", "1⊗0", "1⊗1", "1⊗2"],
                "tensor_product(z2, z3) did not produce the expected labels"
            )

            check_equal(
                fusion_product(tp, 1, 1),
                Dict(1 => 1),
                "tensor_product(z2, z3): vacuum × vacuum was not vacuum"
            )

            # (1⊗0) ⊗ (1⊗0) = 0⊗0
            check_equal(
                fusion_product(tp, 4, 4),
                Dict(1 => 1),
                "tensor_product(z2, z3): (1⊗0) × (1⊗0) was not 0⊗0"
            )
        end

        maybe_testset("reference", "3. reference / Anyonica parity") do
            # TODO
            @test true
        end
    end

    # ============================================================
    # which_permutation / is_equivalent
    # ============================================================

    @testset "which_permutation / is_equivalent" begin
        maybe_testset("basic", "1. basic construction") do
            z3 = zn_fusion_ring(3)
            pz3 = permute(z3, [1,3,2])

            wp = which_permutation(z3, pz3)
            eq = is_equivalent(z3, pz3)

            check_true(wp === missing || wp isa Vector{Int},
                "which_permutation(z3, pz3) returned neither missing nor a permutation vector")
            check_true(eq isa Bool,
                "is_equivalent(z3, pz3) did not return a Bool")
        end

        maybe_testset("intermediate", "2. intermediate correctness") do
            z2 = zn_fusion_ring(2)
            z3 = zn_fusion_ring(3)
            z4 = zn_fusion_ring(4)
            pz3 = permute(z3, [1,3,2])

            wp = which_permutation(z3, pz3)
            check_true(wp !== missing,
                "which_permutation(z3, permute(z3,[1,3,2])) returned missing")
            check_equal(
                permute_mult_tab(multiplication_table(z3), wp),
                multiplication_table(pz3),
                "which_permutation did not return a valid matching permutation"
            )

            check_true(
                is_equivalent(z3, pz3),
                "is_equivalent(z3, permute(z3,[1,3,2])) returned false"
            )

            check_equal(
                which_permutation(z2, z4),
                missing,
                "which_permutation(z2, z4) should have returned missing for different ranks"
            )

            check_false(
                is_equivalent(z2, z4),
                "is_equivalent(z2, z4) should have returned false"
            )
        end

        maybe_testset("reference", "3. reference / Anyonica parity") do
            # TODO
            @test true
        end
    end

end