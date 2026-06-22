using Test
using JSON


#TODO (read to understand how this file works):
#purpose of `test_utils.jl` is to centralize  reusable testing  for  FusionRings.jl test suite. 
# first part defines small wrappers around Julia’s `@test` macros , ex: `check_equal`, `check_approx`, and `check_throws`, 
#failed tests print more informative diagnostic messages before the actual test failure. 
# `check_equal_tensor` is specialized for multiplication tables: instead of printing an entire large 3-dimensional tensor,  reports the expected and actual tensor sizes and the first entry where they differ.
#second part of the file handles the Anyonica reference data. 
# exported Mathematica fixtures all have the same JSON shape, w/ fields `"Input"`, `"Output"`, and `"Info"`,
# file also includes normalization helpers for outputs  such as sub-fusion rings, automorphism lists, tensor-product decompositions, and upper central series data.  
# helpers convert  JSON output and  Julia output into comparable  forms before testing.
#By default, the helper selects a representative subset of cases, while setting `FUSIONRINGS_ORACLE_MODE=all` runs the full Anyonica parity check. 
#how representatives are chosen: always include the first and last cases, and if there are at least 4 cases, also include the cases at 25%, 50%, and 75% through the list.
#1st case
#about 25% through the file
#about 50% through the file
#about 75% through the file
#last case
#ex:
#if there are 2400 cases and max_cases = 12, it picks ~:
#1, 219, 437, 655, 873, 1091, 1200, 1309, 1527, 1745, 1963, 2181, 2400


# ============================================================
# Basic test wrappers
# ============================================================

function check_true(cond::Bool, msg::AbstractString)
    cond || println("FAILED: ", msg)
    @test cond
end

function check_false(cond::Bool, msg::AbstractString)
    cond && println("FAILED: ", msg)
    @test !cond
end

function check_equal(actual, expected, msg::AbstractString)
    if actual != expected
        println("FAILED: ", msg)
        println("  expected: ", repr(expected))
        println("  got     : ", repr(actual))
    end
    @test actual == expected
end

function check_approx(actual, expected, msg::AbstractString; atol = 1e-8, rtol = 1e-8)
    ok = isapprox(actual, expected; atol = atol, rtol = rtol)
    if !ok
        println("FAILED: ", msg)
        println("  expected approximately: ", repr(expected))
        println("  got                   : ", repr(actual))
        println("  absolute difference   : ", repr(abs(actual - expected)))
    end
    @test isapprox(actual, expected; atol = atol, rtol = rtol)
end

function check_throws(f::Function, msg::AbstractString)
    threw = false
    err = nothing

    try
        f()
    catch e
        threw = true
        err = e
    end

    if !threw
        println("FAILED: ", msg)
        println("  expected an exception, but none was thrown")
    end

    @test threw
end

# ============================================================
# Better multiplication-table comparison
# ============================================================

function check_equal_tensor(
    actual::Array{Int,3},
    expected::Array{Int,3},
    msg::AbstractString,
)
    if actual != expected
        println("FAILED: ", msg)
        println("  expected size: ", size(expected))
        println("  actual size  : ", size(actual))

        if size(actual) == size(expected)
            mismatch = findfirst(!=(0), actual .- expected)

            if mismatch !== nothing
                println("  first mismatch at index: ", mismatch)
                println("  expected: ", expected[mismatch])
                println("  got     : ", actual[mismatch])
            else
                println("  arrays differ, but no mismatch was found by subtraction")
            end
        else
            println("  sizes differ, so entrywise comparison was skipped")
        end
    end

    @test actual == expected
end

# ============================================================
# Test-level helpers
# ============================================================

"""
    normalize_test_level(name)

Accept both descriptive names and level-number names.

level1 = basic
level2 = intermediate
level3 = reference
"""
function normalize_test_level(name::AbstractString)
    x = lowercase(strip(name))

    if x in ("1", "level1", "lvl1", "basic")
        return "basic"
    elseif x in ("2", "level2", "lvl2", "intermediate")
        return "intermediate"
    elseif x in ("3", "level3", "lvl3", "reference", "anyonica")
        return "reference"
    elseif x == "all"
        return "all"
    else
        return x
    end
end

# ============================================================
# Anyonica JSON data paths
# ============================================================

"""
Use a normal global instead of `const` so repeated includes/reloads do
not trigger invalid constant redefinition errors.
"""

ANYONICA_TESTDATA_DIR =
    isdefined(@__MODULE__, :ANYONICA_TESTDATA_DIR) ?
    getfield(@__MODULE__, :ANYONICA_TESTDATA_DIR) : joinpath(@__DIR__, "testdata")


"""
    anyonica_data_path(filename)

Find  Anyonica test-data file. Supports  old flat layout:

    test/testdata/zn_tables.json

and new  layout:

    test/testdata/creation/zn_tables.json
    test/testdata/operations/permuted_tabs.json
    test/testdata/properties/fpdim.json
"""
function anyonica_data_path(filename::AbstractString)
    candidates = [
        joinpath(ANYONICA_TESTDATA_DIR, filename),
        joinpath(ANYONICA_TESTDATA_DIR, "creation", filename),
        joinpath(ANYONICA_TESTDATA_DIR, "operations", filename),
        joinpath(ANYONICA_TESTDATA_DIR, "properties", filename),
        joinpath(ANYONICA_TESTDATA_DIR, "formatting", filename),
        joinpath(ANYONICA_TESTDATA_DIR, "import_export", filename),
    ]

    for path in candidates
        if isfile(path)
            return path
        end
    end

    # Return  flat path in  error message because  is the most
    # natural place someone would look first.
    return first(candidates)
end

function has_anyonica_data(filename::AbstractString)
    return isfile(anyonica_data_path(filename))
end

function load_anyonica_data(filename::AbstractString)
    path = anyonica_data_path(filename)
    isfile(path) || error("Missing Anyonica test data file: $path")
    return JSON.parsefile(path)
end

# ============================================================
# JSON conversion helpers
# ============================================================

json_int(x::Integer) = Int(x)
json_int(x::AbstractFloat) = Int(x)
json_int(x::AbstractString) = parse(Int, x)

function json_int_vector(x)
    return Int[json_int(a) for a in x]
end

function json_int_matrix(x)
    m = length(x)
    n = length(x[1])
    A = zeros(Int, m, n)

    for i = 1:m, j = 1:n
        A[i, j] = json_int(x[i][j])
    end

    return A
end

function json_int_3tensor(x)
    r = length(x)
    A = zeros(Int, r, r, r)

    for i = 1:r, j = 1:r, k = 1:r
        A[i, j, k] = json_int(x[i][j][k])
    end

    return A
end

function json_float_vector(x)
    return Float64[Float64(a) for a in x]
end

function json_nested_float_vector(x)
    return [json_float_vector(row) for row in x]
end

# ============================================================
# Detecting JSON object shapes
# ============================================================

function looks_like_int_vector(x)
    x isa AbstractVector || return false
    isempty(x) && return true
    return all(a -> a isa Integer || a isa AbstractFloat || a isa AbstractString, x)
end

function looks_like_matrix(x)
    x isa AbstractVector || return false
    isempty(x) && return false
    x[1] isa AbstractVector || return false
    isempty(x[1]) && return false
    return !(x[1][1] isa AbstractVector)
end

function looks_like_3tensor(x)
    x isa AbstractVector || return false
    isempty(x) && return false
    x[1] isa AbstractVector || return false
    isempty(x[1]) && return false
    x[1][1] isa AbstractVector || return false
    return true
end

"""
    looks_like_formal_code(x)

Formal codes look like short integer vectors `[rank, multiplicity, nnsd, index]`.
"""
function looks_like_formal_code(x)
    return looks_like_int_vector(x)
end

# ============================================================
# Ring construction from Anyonica references
# ============================================================

function ring_from_anyonica_mt(x)
    return fusion_ring(json_int_3tensor(x); skip_check = true)
end

function code_from_json(x)
    return json_int_vector(x)
end

function ring_from_anyonica_code(x)
    code = code_from_json(x)
    return from_anyonwiki_code(code)
end

"""
    ring_from_anyonica_ref(x)

Accept either:
- a raw multiplication table, or
- a formal AnyonWiki code.

This supports both old and new fixture formats.
"""
function ring_from_anyonica_ref(x)
    if looks_like_3tensor(x)
        return ring_from_anyonica_mt(x)
    elseif looks_like_formal_code(x)
        return ring_from_anyonica_code(x)
    else
        error(
            "Could not interpret Anyonica ring reference of type $(typeof(x)): $(repr(x))",
        )
    end
end

# Alias with a short name for tests.
ring_from_ref(x) = ring_from_anyonica_ref(x)

# ============================================================
# Multiplication-table oracle comparison
# ============================================================

function mt_from_anyonica_ref(x)
    if looks_like_3tensor(x)
        return json_int_3tensor(x)
    elseif looks_like_formal_code(x)
        return multiplication_table(ring_from_anyonica_code(x))
    else
        error("Could not interpret multiplication-table reference: $(repr(x))")
    end
end

function check_mt_equal(actual_mt::Array{Int,3}, expected_json, msg::AbstractString)
    expected = json_int_3tensor(expected_json)
    check_equal_tensor(actual_mt, expected, msg)
end

function check_mt_equal(actual_ring::FusionRing, expected_json, msg::AbstractString)
    check_mt_equal(multiplication_table(actual_ring), expected_json, msg)
end

# ============================================================
# Permutation helpers
# ============================================================

function first_permutation_result(x)
    x === nothing && return nothing
    x === missing && return nothing
    x isa Vector{Int} && return x
    x isa Vector{Vector{Int}} && return isempty(x) ? nothing : first(x)
    return x
end

function nontrivial_permutations(perms::Vector{Vector{Int}}, r::Int)
    id = collect(1:r)
    return sort([p for p in perms if p != id])
end

function sorted_json_permutations(xs)
    return sort([json_int_vector(p) for p in xs])
end

# ============================================================
# Subring normalization helpers
# ============================================================

"""
    normalize_subrings_from_anyonica(xs)



old:
    [els, multiplication_table]

new:
    [els, formal_code]
"""
function normalize_subrings_from_anyonica(xs)
    out = Tuple{Vector{Int},Array{Int,3}}[]

    for pair in xs
        els = json_int_vector(pair[1])
        sub_ref = pair[2]
        sub_mt = mt_from_anyonica_ref(sub_ref)
        push!(out, (els, sub_mt))
    end

    return sort(out; by = x -> x[1])
end

function normalize_subring_sets_from_anyonica(xs)
    return sort([json_int_vector(pair[1]) for pair in xs])
end

function normalize_subrings_from_julia(xs)
    out = Tuple{Vector{Int},Array{Int,3}}[]

    for item in xs
        if item isa Tuple
            els = collect(Int, item[1])
            subring = item[2]
            push!(out, (els, multiplication_table(subring)))
        elseif item isa Dict
            inj =
                haskey(item, "injection") ? item["injection"] :
                haskey(item, :injection) ? item[:injection] :
                error("Subring Dict did not contain an injection key.")

            sub =
                haskey(item, "fusion_ring") ? item["fusion_ring"] :
                haskey(item, :fusion_ring) ? item[:fusion_ring] :
                error("Subring Dict did not contain a fusion_ring key.")

            push!(out, (collect(Int, inj), multiplication_table(sub)))
        else
            error("Unsupported sub_fusion_rings output item type: $(typeof(item))")
        end
    end

    return sort(out; by = x -> x[1])
end

# Code/pair normalization helpers

function normalize_code(x)
    return json_int_vector(x)
end

function normalize_code_list(xs)
    return sort([normalize_code(x) for x in xs])
end

function normalize_decomposition_codes(xs)
    # xs  usually  list of decompositions, each decomposition being  list of factor codes.
    return sort([sort([normalize_code(code) for code in dec]) for dec in xs])
end

function normalize_pair_code_output(x)
    # For outputs like [els, code].
    return (json_int_vector(x[1]), normalize_code(x[2]))
end

function normalize_series_code_output(xs)
    # For outputs like [[els1, code1], [els2, code2], ...].
    return [(json_int_vector(pair[1]), normalize_code(pair[2])) for pair in xs]
end


"""
Anyonica numerical outputs may be  numbers, strings, or nested lists.
 converts simple numeric-ish values to Float64.
"""
function json_float(x::Integer)
    return Float64(x)
end

function json_float(x::AbstractFloat)
    return Float64(x)
end

function json_float(x::AbstractString)
    return parse(Float64, x)
end

function numeric_vector_from_json(xs)
    return Float64[json_float(x) for x in xs]
end

function check_numeric_vector_approx(
    actual,
    expected_json,
    msg::AbstractString;
    atol = 1e-8,
    rtol = 1e-8,
)
    actual_vec = Float64[Float64(real(x)) for x in actual]
    expected_vec = numeric_vector_from_json(expected_json)

    if length(actual_vec) != length(expected_vec)
        println("FAILED: ", msg)
        println("  expected length: ", length(expected_vec))
        println("  actual length  : ", length(actual_vec))
    elseif !isapprox(actual_vec, expected_vec; atol = atol, rtol = rtol)
        println("FAILED: ", msg)
        println("  expected approximately: ", repr(expected_vec))
        println("  got                   : ", repr(actual_vec))

        diffs = abs.(actual_vec .- expected_vec)
        worst = argmax(diffs)
        println("  worst mismatch index  : ", worst)
        println("  expected there        : ", expected_vec[worst])
        println("  got there             : ", actual_vec[worst])
        println("  absolute difference   : ", diffs[worst])
    end

    @test length(actual_vec) == length(expected_vec)
    @test isapprox(actual_vec, expected_vec; atol = atol, rtol = rtol)
end

function check_numeric_approx(
    actual,
    expected_json,
    msg::AbstractString;
    atol = 1e-8,
    rtol = 1e-8,
)
    actual_num = Float64(real(actual))
    expected_num = json_float(expected_json)

    if !isapprox(actual_num, expected_num; atol = atol, rtol = rtol)
        println("FAILED: ", msg)
        println("  expected approximately: ", repr(expected_num))
        println("  got                   : ", repr(actual_num))
        println("  absolute difference   : ", repr(abs(actual_num - expected_num)))
    end

    @test isapprox(actual_num, expected_num; atol = atol, rtol = rtol)
end

# Oracle case selection

"""
    oracle_mode()

Controls how much JSON oracle data to test.

Default:
- If FUSIONRINGS_ORACLE_MODE is unset, use `representative`.
- Set FUSIONRINGS_ORACLE_MODE=all for full Anyonica parity.
"""
function oracle_mode()
    mode = lowercase(get(ENV, "FUSIONRINGS_ORACLE_MODE", "representative"))

    if mode in ("all", "full")
        return :all
    elseif mode == "first"
        return :first
    elseif mode == "spread"
        return :spread
    elseif mode in ("representative", "rep")
        return :representative
    elseif mode in ("small", "small_ranks")
        return :small_ranks
    else
        error("Unknown FUSIONRINGS_ORACLE_MODE=$mode")
    end
end

function _rank_from_oracle_input(x)
    if looks_like_formal_code(x)
        return json_int(x[1])
    elseif looks_like_3tensor(x)
        return length(x)
    else
        return nothing
    end
end

function oracle_case_indices(data; mode = oracle_mode(), max_cases = 12, max_rank = 6)
    n = length(data["Input"])

    n == 0 && return Int[]

    if mode == :all
        return collect(1:n)
    elseif mode == :first
        return collect(1:min(max_cases, n))
    elseif mode == :spread
        if n <= max_cases
            return collect(1:n)
        end
        return sort(unique(round.(Int, range(1, n; length = max_cases))))
    elseif mode == :representative
        if n <= max_cases
            return collect(1:n)
        end

        base = Int[1, max(1, cld(n, 4)), max(1, cld(n, 2)), max(1, cld(3n, 4)), n]

        spread = round.(Int, range(1, n; length = max_cases))
        return sort(unique(vcat(base, collect(spread))))
    elseif mode == :small_ranks
        inds = Int[]

        for i = 1:n
            r = _rank_from_oracle_input(data["Input"][i])
            if r !== nothing && r <= max_rank
                push!(inds, i)
            end
        end

        return inds
    else
        error("Unknown oracle selection mode: $mode")
    end
end

function oracle_case_indices(
    filename::AbstractString,
    data;
    mode = oracle_mode(),
    max_cases = 12,
    max_rank = 6,
)
    # Creation  are tiny, so testing all of them is cheap and useful.
    if filename in (
        "zn_tables.json",
        "grouptables.json",
        "HI_tables.json",
        "TY_tables.json",
        "psu2k_tables.json",
        "su2k_tables.json",
        "son2_tables.json",
    )
        return collect(1:length(data["Input"]))
    end

    return oracle_case_indices(
        data;
        mode = mode,
        max_cases = max_cases,
        max_rank = max_rank,
    )
end
