@testset "Import / export" begin
    #TODO: changes to import/export tests:
    #Added reference tests for mttojs and ring_to_dict
    #Used multiplication-table outputs from:
    #tensor_product_tables.json
    #to build FusionRing objects, export them to temporary JSON files, import them back, and compare multiplication tables.
    #For full multiplication tables, replaced ordinary check_equal with:
    #check_equal_tensor
    #Used creation/zn_tables.json
    #to verify that Anyonica-style nested multiplication tables can be decoded with:
    #mtfromjs(Dict("mult_tab" => ...))
    #Used rings from:fpdim.json
    #to check:
    #from_qqb_id(qqb_id(fpdim(r))) == fpdim(r)
    #Uses Anyonica multiplication tables from:
    #zn_tables.json
    #to build rings and check:

    #mttojs(r)
    #ring_to_dict(r)
    #mtfromjs(ring_to_dict(r))

  # ============================================================
  # QQBar ID helpers
  # ============================================================

    @testset "qqb_id / from_qqb_id" begin
        maybe_testset("basic", "1. basic construction") do
            r = zn_fusion_ring(2)
            d = fpdim(r)
            s = qqb_id(d)

            check_true(
                s isa AbstractString,
                "qqb_id(fpdim(zn_fusion_ring(2))) did not return a string",
            )
            check_false(
                isempty(s),
                "qqb_id(fpdim(zn_fusion_ring(2))) returned an empty string",
            )
        end

    maybe_testset("intermediate", "2. intermediate correctness") do
      r = zn_fusion_ring(2)
      d = fpdim(r)
      s = qqb_id(d)
      d2 = from_qqb_id(s)

            check_equal(
                d2,
                d,
                "from_qqb_id(qqb_id(x)) did not recover x for fpdim(zn_fusion_ring(2))",
            )

      ds = fpdims(r)
      ids = qqb_id(ds)
      ds2 = from_qqb_id(Any[ids...])

            check_equal(
                ds2,
                ds,
                "from_qqb_id.(qqb_id.(fpdims(r))) did not recover the FP dimensions",
            )
        end

        maybe_testset("reference", "3. reference / Anyonica parity") do
            data = load_anyonica_data("fpdim.json")

            for idx in oracle_case_indices("fpdim.json", data; max_cases = 8)
                r = ring_from_anyonica_code(data["Input"][idx])
                d = fpdim(r)
                s = qqb_id(d)
                d2 = from_qqb_id(s)

                check_equal(
                    d2,
                    d,
                    "from_qqb_id(qqb_id(fpdim(r))) failed for Anyonica fpdim.json case $idx",
                )
            end
        end
    end

  # ============================================================
  # Low-level JSON decoding helpers
  # ============================================================

    @testset "low-level JSON decoding helpers" begin
        maybe_testset("basic", "1. basic construction") do
            js_fc = Dict("formal_code" => [2, 1, 0, 7])
            js_mt = Dict("mult_tab" => [[[1, 0], [0, 1]], [[0, 1], [1, 0]]])
            js_bc = Dict("barcode" => "123456789")
            js_tpd = Dict("tensor_product_decompositions" => Any[])
            js_nch = Dict("numeric_characters" => [[[1.0, 0.0]]])
            js_nfpds =
                Dict("numeric_frobenius_perron_dimensions" => [[1.0, 0.0], [2.0, 0.0]])
            js_nfpd = Dict("numeric_frobenius_perron_dimension" => [3.0, 0.0])

            check_true(
                fcfromjs(js_fc) isa Vector,
                "fcfromjs did not return a vector on a basic formal_code input",
            )
            check_true(
                mtfromjs(js_mt) isa Array{Int,3},
                "mtfromjs did not return an Int 3-tensor on a basic mult_tab input",
            )
            check_true(
                bcfromjs(js_bc) isa ZZRingElem,
                "bcfromjs did not return a ZZRingElem on a basic barcode input",
            )
            check_true(
                tpdfromjs(js_tpd) isa Vector,
                "tpdfromjs did not return a vector on an empty decomposition input",
            )
            check_true(
                nchfromjs(js_nch) isa Matrix{ComplexF64},
                "nchfromjs did not return a ComplexF64 matrix on a basic input",
            )
            check_true(
                nfpdsfromjs(js_nfpds) isa Vector{ComplexF64},
                "nfpdsfromjs did not return a ComplexF64 vector on a basic input",
            )
            check_true(
                nfpdfromjs(js_nfpd) isa ComplexF64,
                "nfpdfromjs did not return a ComplexF64 on a basic input",
            )
        end

    maybe_testset("intermediate", "2. intermediate correctness") do
      js_fc1 = Dict("formal_code" => [2, 1, 0, 7])
      js_fc2 = Dict("anyonwiki_code" => [3, 1, 2, 9])
      js_fc3 = Dict("formal_code" => Any[])

            check_equal(
                fcfromjs(js_fc1),
                [2, 1, 0, 7],
                "fcfromjs did not decode formal_code correctly",
            )
            check_equal(
                fcfromjs(js_fc2),
                [3, 1, 2, 9],
                "fcfromjs did not decode legacy anyonwiki_code correctly",
            )
            check_equal(
                fcfromjs(js_fc3),
                missing,
                "fcfromjs did not return missing on an empty code",
            )

            js_mt = Dict("mult_tab" => [[[1, 0], [0, 1]], [[0, 1], [1, 0]]])
            expected_mt = zeros(Int, 2, 2, 2)
            expected_mt[1, 1, 1] = 1
            expected_mt[1, 2, 2] = 1
            expected_mt[2, 1, 2] = 1
            expected_mt[2, 2, 1] = 1

            check_equal_tensor(
                mtfromjs(js_mt),
                expected_mt,
                "mtfromjs did not reconstruct the expected Z2 multiplication table",
            )

            js_bc = Dict("barcode" => "123456789")
            check_equal(
                string(bcfromjs(js_bc)),
                "123456789",
                "bcfromjs did not reconstruct the expected barcode value",
            )

            js_tpd1 = Dict("tensor_product_decompositions" => Any[])
            check_equal(
                tpdfromjs(js_tpd1),
                Any[],
                "tpdfromjs did not return [] for an empty decomposition list",
            )

            js_tpd2 = Dict(
                "tensor_product_decompositions" =>
                    [[[2, 1, 0, 1], [3, 1, 0, 1]], [[4, 1, 0, 1]]],
            )
            check_equal(
                tpdfromjs(js_tpd2),
                [[[2, 1, 0, 1], [3, 1, 0, 1]], [[4, 1, 0, 1]]],
                "tpdfromjs did not decode direct decomposition data correctly",
            )

            js_nch1 = Dict("numeric_characters" => nothing)
            check_equal(
                nchfromjs(js_nch1),
                missing,
                "nchfromjs did not return missing when numeric_characters was nothing",
            )

            js_nch2 = Dict(
                "numeric_characters" =>
                    [[[1.0, 0.0], [0.0, 1.0]], [[0.0, -1.0], [2.0, 0.5]]],
            )
            expected_nch2 = ComplexF64[
                1.0 + 0.0im 0.0 + 1.0im
                0.0 - 1.0im 2.0 + 0.5im
            ]
            check_equal(
                nchfromjs(js_nch2),
                expected_nch2,
                "nchfromjs did not decode a 2×2 complex matrix correctly",
            )

            js_nfpds =
                Dict("numeric_frobenius_perron_dimensions" => [[1.0, 0.0], [2.5, -1.0]])
            check_equal(
                nfpdsfromjs(js_nfpds),
                ComplexF64[1.0+0.0im, 2.5-1.0im],
                "nfpdsfromjs did not decode numeric FP dimensions correctly",
            )

            js_nfpd = Dict("numeric_frobenius_perron_dimension" => [3.0, -2.0])
            check_equal(
                nfpdfromjs(js_nfpd),
                3.0 - 2.0im,
                "nfpdfromjs did not decode numeric FP dimension correctly",
            )

            js_names = Dict("names" => ["A", "B"])
            js_texnames = Dict("texnames" => ["A", "B"])
            check_equal(
                nfromjs(js_names),
                ["A", "B"],
                "nfromjs did not decode names correctly",
            )
            check_equal(
                tnfromjs(js_texnames),
                ["A", "B"],
                "tnfromjs did not decode texnames correctly",
            )
        end

        maybe_testset("reference", "3. reference / Anyonica parity") do
            # Creation fixtures give raw multiplication-table outputs.
            # This tests that mtfromjs can decode Anyonica-style nested mult_tab data.
            data = load_anyonica_data("zn_tables.json")

            for idx in oracle_case_indices("zn_tables.json", data)
                expected = json_int_3tensor(data["Output"][idx])
                js = Dict("mult_tab" => data["Output"][idx])

                check_equal_tensor(
                    mtfromjs(js),
                    expected,
                    "mtfromjs failed to decode Anyonica zn_tables.json output case $idx",
                )
            end

            # Properties fixtures use formal codes.  This tests both supported keys:
            # formal_code and legacy anyonwiki_code.
            props = load_anyonica_data("is_group_ring.json")

            for idx in oracle_case_indices("is_group_ring.json", props; max_cases = 8)
                code = json_int_vector(props["Input"][idx])

                check_equal(
                    fcfromjs(Dict("formal_code" => code)),
                    code,
                    "fcfromjs failed on formal_code for Anyonica properties case $idx",
                )
                check_equal(
                    fcfromjs(Dict("anyonwiki_code" => code)),
                    code,
                    "fcfromjs failed on anyonwiki_code for Anyonica properties case $idx",
                )
            end
        end
    end

  # ============================================================
  # Misc import-side helpers
  # ============================================================

    @testset "misc import-side helpers" begin
        maybe_testset("basic", "1. basic construction") do
            check_true(
                vec_to_cflt([1.0, 2.0]) isa ComplexF64,
                "vec_to_cflt([re,im]) did not return a ComplexF64",
            )

            js_false = Dict("categorifiable" => false)
            js_null = Dict("categorifiable" => nothing)
            js_true = Dict("categorifiable" => true)

            check_true(
                cfromjs(js_false) isa Bool,
                "cfromjs(false case) did not return a Bool",
            )
            check_true(
                cfromjs(js_true) isa Bool,
                "cfromjs(true case) did not return a Bool",
            )
        end

        maybe_testset("intermediate", "2. intermediate correctness") do
            check_equal(
                vec_to_cflt([1.0, 2.0]),
                1.0 + 2.0im,
                "vec_to_cflt([1.0,2.0]) did not decode correctly",
            )

            check_equal(
                cfromjs(Dict("categorifiable" => false)),
                false,
                "cfromjs did not decode categorifiable=false correctly",
            )
            check_equal(
                cfromjs(Dict("categorifiable" => nothing)),
                missing,
                "cfromjs did not decode categorifiable=nothing correctly",
            )
            check_equal(
                cfromjs(Dict("categorifiable" => true)),
                true,
                "cfromjs did not decode categorifiable=true correctly",
            )

            check_equal(
                ctsfromjs(Dict("categorifications" => nothing)),
                missing,
                "ctsfromjs did not return missing when categorifications was nothing",
            )
            check_equal(
                ctsfromjs(Dict("categorifications" => Any[])),
                Vector{Int64}[],
                "ctsfromjs did not return an empty vector-of-vectors on empty input",
            )
            check_equal(
                ctsfromjs(Dict("categorifications" => [[1, 2, 3, 4], [5, 6, 7, 8]])),
                [[1, 2, 3, 4], [5, 6, 7, 8]],
                "ctsfromjs did not decode categorification codes correctly",
            )

            check_equal(
                ncrfromjs(Dict("non_cat_reasons" => Dict("Fusion" => "x"))),
                Dict("Fusion" => "x"),
                "ncrfromjs did not return explicit non_cat_reasons when present",
            )

            default_ncr = ncrfromjs(Dict())
            check_true(
                haskey(default_ncr, "Fusion"),
                "ncrfromjs default dictionary did not contain key \"Fusion\"",
            )
            check_true(
                haskey(default_ncr, "Modular"),
                "ncrfromjs default dictionary did not contain key \"Modular\"",
            )
        end

        maybe_testset("reference", "3. reference / Anyonica parity") do
            data = load_anyonica_data("fpdims.json")

            for idx in oracle_case_indices("fpdims.json", data; max_cases = 8)
                expected = numeric_vector_from_json(data["Output"][idx])
                encoded = [[x, 0.0] for x in expected]
                decoded =
                    nfpdsfromjs(Dict("numeric_frobenius_perron_dimensions" => encoded))

                check_equal(
                    decoded,
                    ComplexF64[x + 0.0im for x in expected],
                    "nfpdsfromjs failed to decode Anyonica-style FP dimensions case $idx",
                )
            end
        end
    end

  # ============================================================
  # Export-side helpers
  # ============================================================

  @testset "export-side helpers" begin
    maybe_testset("basic", "1. basic construction") do
      r = zn_fusion_ring(2)

            check_true(
                mttojs(r) isa Vector,
                "mttojs(zn_fusion_ring(2)) did not return a nested vector",
            )
            check_true(
                fusion_ring_string(r) isa AbstractString,
                "fusion_ring_string(zn_fusion_ring(2)) did not return a string",
            )
            check_true(
                fusion_ring_file_name(r) isa AbstractString,
                "fusion_ring_file_name(zn_fusion_ring(2)) did not return a string",
            )
            check_true(
                ring_to_dict(r) isa Dict,
                "ring_to_dict(zn_fusion_ring(2)) did not return a Dict",
            )
        end

    maybe_testset("intermediate", "2. intermediate correctness") do
      r = zn_fusion_ring(2)

            expected_mt = [[[1, 0], [0, 1]], [[0, 1], [1, 0]]]
            check_equal(
                mttojs(r),
                expected_mt,
                "mttojs(zn_fusion_ring(2)) did not produce the expected nested representation",
            )

            check_equal(
                missing_to_nothing(missing),
                nothing,
                "missing_to_nothing(missing) was not nothing",
            )
            check_equal(missing_to_nothing(5), 5, "missing_to_nothing(5) did not return 5")

            check_equal(
                reim(1.5 + 2.0im),
                [1.5, 2.0],
                "reim(::ComplexF64) did not split real/imag parts correctly",
            )
            check_equal(reim(3.0), [3.0, 0.0], "reim(::Float64) did not produce [x,0.0]")

            mat = ComplexF64[1.0+0.0im 2.0-1.0im; 0.0+3.0im -1.5+0.5im]
            check_equal(
                reim(mat),
                [[[1.0, 0.0], [2.0, -1.0]], [[0.0, 3.0], [-1.5, 0.5]]],
                "reim(::Matrix{ComplexF64}) did not split entries correctly",
            )

            vv = [[1.0+0.0im, 2.0-1.0im], [0.0+3.0im]]
            check_equal(
                reim(vv),
                [[[1.0, 0.0], [2.0, -1.0]], [[0.0, 3.0]]],
                "reim(::Vector{Vector{ComplexF64}}) did not split entries correctly",
            )

            check_true(
                endswith(fusion_ring_file_name(r), ".json"),
                "fusion_ring_file_name(zn_fusion_ring(2)) did not end in .json",
            )

            d = ring_to_dict(r)
            check_true(
                haskey(d, "mult_tab"),
                "ring_to_dict(zn_fusion_ring(2)) did not contain key \"mult_tab\"",
            )
            check_true(
                haskey(d, "anyonwiki_code"),
                "ring_to_dict(zn_fusion_ring(2)) did not contain key \"anyonwiki_code\"",
            )
            check_true(
                haskey(d, "frobenius_perron_dimension"),
                "ring_to_dict(zn_fusion_ring(2)) did not contain key \"frobenius_perron_dimension\"",
            )
        end

        maybe_testset("reference", "3. reference / Anyonica parity") do
            data = load_anyonica_data("zn_tables.json")

            for idx in oracle_case_indices("zn_tables.json", data)
                r = fusion_ring(json_int_3tensor(data["Output"][idx]); skip_check = true)

                check_equal(
                    mttojs(r),
                    data["Output"][idx],
                    "mttojs did not reproduce Anyonica zn_tables.json output case $idx",
                )

                d = ring_to_dict(r)
                check_true(
                    haskey(d, "mult_tab"),
                    "ring_to_dict output did not include mult_tab for Anyonica zn_tables.json case $idx",
                )
                check_equal_tensor(
                    mtfromjs(d),
                    multiplication_table(r),
                    "mtfromjs(ring_to_dict(r)) failed for Anyonica zn_tables.json case $idx",
                )
            end
        end
    end

  # ============================================================
  # End-to-end import/export
  # ============================================================

  @testset "end-to-end import/export" begin
    maybe_testset("basic", "1. basic construction") do
      mktempdir() do dir
        r = zn_fusion_ring(2)
        file = joinpath(dir, "ring.json")

                export_ring(file, r)
                check_true(isfile(file), "export_ring did not create the output file")

                r2 = import_ring(file)
                check_true(r2 isa FusionRing, "import_ring did not return a FusionRing")
            end
        end

    maybe_testset("intermediate", "2. intermediate correctness") do
      mktempdir() do dir
        r = zn_fusion_ring(2)
        file = joinpath(dir, "ring.json")

        export_ring(file, r)
        r2 = import_ring(file)

                check_equal_tensor(
                    multiplication_table(r2),
                    multiplication_table(r),
                    "import_ring(export_ring(r)) did not preserve the multiplication table",
                )
                check_equal(
                    rank(r2),
                    rank(r),
                    "import_ring(export_ring(r)) did not preserve the rank",
                )
                check_equal(
                    labels(r2),
                    labels(r),
                    "import_ring(export_ring(r)) did not preserve labels",
                )
            end
        end

        maybe_testset("reference", "3. reference / Anyonica parity") do
            data = load_anyonica_data("tensor_product_tables.json")

            for idx in
                oracle_case_indices("tensor_product_tables.json", data; max_cases = 6)
                r = fusion_ring(json_int_3tensor(data["Output"][idx]); skip_check = true)

                mktempdir() do dir
                    file = joinpath(dir, "anyonica_tensor_product_case_$idx.json")

                    export_ring(file, r)
                    r2 = import_ring(file)

                    check_equal_tensor(
                        multiplication_table(r2),
                        multiplication_table(r),
                        "import_ring(export_ring(r)) failed for Anyonica tensor_product_tables.json output case $idx",
                    )
                    check_equal(
                        rank(r2),
                        rank(r),
                        "import_ring(export_ring(r)) did not preserve rank for Anyonica tensor_product_tables.json output case $idx",
                    )
                end
            end
        end
    end

  # ============================================================
  # Multi-ring export/import
  # ============================================================

  @testset "multi-ring export/import" begin
    maybe_testset("basic", "1. basic construction") do
      rs = [zn_fusion_ring(2), zn_fusion_ring(3)]
      d = rings_to_dict(rs)

            check_true(d isa Dict, "rings_to_dict did not return a Dict")
            check_true(
                haskey(d, "data"),
                "rings_to_dict output did not contain key \"data\"",
            )

            mktempdir() do dir
                file = joinpath(dir, "rings.json")
                export_rings(file, rs)
                check_true(isfile(file), "export_rings did not create the output file")
            end
        end

    maybe_testset("intermediate", "2. intermediate correctness") do
      rs = [zn_fusion_ring(2), zn_fusion_ring(3)]
      d = rings_to_dict(rs)

            check_true(
                haskey(d, "info"),
                "rings_to_dict output did not contain key \"info\"",
            )

            mktempdir() do dir
                file = joinpath(dir, "rings.json")
                export_rings(file, rs)

                loaded = import_rings(file)
                check_equal(
                    length(loaded),
                    length(rs),
                    "import_rings(export_rings(rs)) did not preserve the number of rings",
                )

                for i in eachindex(rs)
                    check_equal_tensor(
                        multiplication_table(loaded[i]),
                        multiplication_table(rs[i]),
                        "import_rings(export_rings(rs)) did not preserve multiplication table for ring $i",
                    )
                end
            end
        end

        maybe_testset("reference", "3. reference / Anyonica parity") do
            data = load_anyonica_data("is_group_ring.json")
            inds = oracle_case_indices("is_group_ring.json", data; max_cases = 5)

            rs = [ring_from_anyonica_code(data["Input"][idx]) for idx in inds]

            mktempdir() do dir
                file = joinpath(dir, "anyonica_rings.json")

                export_rings(file, rs)
                loaded = import_rings(file)

                check_equal(
                    length(loaded),
                    length(rs),
                    "import_rings(export_rings(rs)) did not preserve number of Anyonica reference rings",
                )

                for i in eachindex(rs)
                    check_equal_tensor(
                        multiplication_table(loaded[i]),
                        multiplication_table(rs[i]),
                        "import_rings(export_rings(rs)) did not preserve Anyonica reference ring $i",
                    )
                end
            end
        end
    end
end
