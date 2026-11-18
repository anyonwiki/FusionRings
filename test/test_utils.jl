using Test
using JSON

# Shared assertions, fixture decoding, and deterministic reference sampling.
# A reference match is parity evidence; it is not independent theorem validation.

# ============================================================
# Assertions
# ============================================================

function check_true(cond, msg::AbstractString)
  ok = cond === true
  ok || println("FAILED: ", msg, "\n  expected true; got ", repr(cond))
  @test ok
  return nothing
end

function check_false(cond, msg::AbstractString)
  ok = cond === false
  ok || println("FAILED: ", msg, "\n  expected false; got ", repr(cond))
  @test ok
  return nothing
end

function check_equal(actual, expected, msg::AbstractString)
  # isequal handles missing metadata without producing a non-Boolean result.
  ok = isequal(actual, expected)
  if !ok
    println("FAILED: ", msg)
    println("  expected: ", repr(expected))
    println("  got     : ", repr(actual))
  end
  @test ok
  return nothing
end

function check_approx(actual, expected, msg::AbstractString; atol = 1e-8, rtol = 1e-8)
  if actual isa AbstractArray && expected isa AbstractArray
    same_axes = axes(actual) == axes(expected)
    if !same_axes
      println("FAILED: ", msg)
      println("  expected axes: ", axes(expected))
      println("  actual axes  : ", axes(actual))
      @test same_axes
      return nothing
    end
  end

  ok = isapprox(actual, expected; atol = atol, rtol = rtol)
  if !ok
    println("FAILED: ", msg)
    println("  expected approximately: ", repr(expected))
    println("  got                   : ", repr(actual))
    println("  atol / rtol           : ", atol, " / ", rtol)

    if actual isa Number && expected isa Number
      println("  absolute difference   : ", abs(actual - expected))
    elseif actual isa AbstractArray && expected isa AbstractArray
      if !isempty(actual)
        println("  max absolute difference: ", maximum(abs.(actual .- expected)))
      end
    end
  end

  @test ok
  return nothing
end

"""
    check_throws(f, ExpectedException, description; message_contains=nothing)

Check the exception type and, optionally, a stable substring of its message.
The expression runs once. InterruptException is always rethrown.
Also supports do-block syntax: check_throws(ArgumentError, description) do ... end.
"""
function check_throws(
  f::Function,
  expected::Type{<:Exception},
  msg::AbstractString;
  message_contains::Union{Nothing,AbstractString} = nothing,
)
  caught = nothing
  did_throw = false

  try
    f()
  catch err
    err isa InterruptException && rethrow()

    # Broad legacy checks must not pass because a name was misspelled.
    # An intentional UndefVarError test can request that exact type explicitly.
    expected === Exception && err isa UndefVarError && rethrow()

    caught = err
    did_throw = true
  end

  type_ok = did_throw && isa(caught, expected)
  if !type_ok
    println("FAILED: ", msg)
    println("  expected exception: ", expected)
    if did_throw
      println("  actual exception  : ", typeof(caught))
      println("  message           : ", sprint(showerror, caught))
    else
      println("  no exception was thrown")
    end
  end
  @test type_ok

  if type_ok && message_contains !== nothing
    text = sprint(showerror, caught)
    message_ok = occursin(message_contains, text)
    if !message_ok
      println("FAILED: ", msg)
      println("  expected message to contain: ", repr(message_contains))
      println("  actual message            : ", text)
    end
    @test message_ok
  end

  return nothing
end

# Temporary compatibility while the remaining test files are revised.
# This form still accepts different exception types, including MethodError:
# the current Float64-constructor rejection test legitimately expects dispatch
# failure. Each call site must choose its intended type when its file is revised.
function check_throws(f::Function, msg::AbstractString)
  @warn "Untyped check_throws call: specify the expected exception type when revising this test." maxlog = 1
  return check_throws(f, Exception, msg)
end

function check_equal_tensor(
  actual::AbstractArray{<:Integer,3},
  expected::AbstractArray{<:Integer,3},
  msg::AbstractString,
)
  same_axes = axes(actual) == axes(expected)
  ok = same_axes && actual == expected

  if !ok
    println("FAILED: ", msg)
    println("  expected size: ", size(expected))
    println("  actual size  : ", size(actual))

    if same_axes
      # Compare entries directly: no subtraction, overflow, or difference tensor.
      for index in CartesianIndices(actual)
        if actual[index] != expected[index]
          println("  first mismatch at: ", index)
          println("  expected         : ", expected[index])
          println("  got              : ", actual[index])
          break
        end
      end
    else
      println("  expected axes: ", axes(expected))
      println("  actual axes  : ", axes(actual))
    end
  end

  @test ok
  return nothing
end

# ============================================================
# Test-level aliases (runtests.jl validates the normalized value)
# ============================================================

function normalize_test_level(name::AbstractString)
  level = lowercase(strip(name))
  level in ("1", "level1", "lvl1", "basic") && return "basic"
  level in ("2", "level2", "lvl2", "intermediate") && return "intermediate"
  level in ("3", "level3", "lvl3", "reference", "anyonica") && return "reference"
  return level
end

# ============================================================
# Fixture lookup and envelope validation
# ============================================================

ANYONICA_TESTDATA_DIR = if isdefined(@__MODULE__, :ANYONICA_TESTDATA_DIR)
  getfield(@__MODULE__, :ANYONICA_TESTDATA_DIR)
else
  joinpath(@__DIR__, "testdata")
end

function anyonica_data_path(filename::AbstractString)
  candidates = unique(abspath.([
    joinpath(ANYONICA_TESTDATA_DIR, filename),
    joinpath(ANYONICA_TESTDATA_DIR, "creation", filename),
    joinpath(ANYONICA_TESTDATA_DIR, "operations", filename),
    joinpath(ANYONICA_TESTDATA_DIR, "properties", filename),
    joinpath(ANYONICA_TESTDATA_DIR, "formatting", filename),
    joinpath(ANYONICA_TESTDATA_DIR, "import_export", filename),
  ]))
  matches = filter(isfile, candidates)

  length(matches) > 1 && throw(ArgumentError(
    "Ambiguous Anyonica fixture $(repr(filename)); matching files: " *
    join(matches, ", "),
  ))

  return isempty(matches) ? first(candidates) : only(matches)
end

has_anyonica_data(filename::AbstractString) = isfile(anyonica_data_path(filename))

function validate_anyonica_data(data, filename::AbstractString)
  data isa AbstractDict || throw(ArgumentError(
    "Invalid fixture $filename: expected a JSON object.",
  ))

  for key in ("Input", "Output", "Info")
    haskey(data, key) || throw(ArgumentError(
      "Invalid fixture $filename: missing key $(repr(key)).",
    ))
  end

  data["Input"] isa AbstractVector || throw(ArgumentError(
    "Invalid fixture $filename: Input must be an array.",
  ))
  data["Output"] isa AbstractVector || throw(ArgumentError(
    "Invalid fixture $filename: Output must be an array.",
  ))
  data["Info"] isa AbstractString || throw(ArgumentError(
    "Invalid fixture $filename: Info must be a string.",
  ))

  length(data["Input"]) == length(data["Output"]) || throw(ArgumentError(
    "Invalid fixture $filename: Input and Output lengths differ.",
  ))
  isempty(data["Input"]) && throw(ArgumentError(
    "Invalid fixture $filename: no reference cases are present.",
  ))

  # Guard against the confirmed exporter bug: this file currently contains
  # universal-grading outputs instead of four criterion Booleans.
  if basename(filename) == "categorifiability.json"
    for (index, output) in enumerate(data["Output"])
      # Preserve explicit unavailable results; the domain test must report skips.
      output === nothing && continue
      valid = output isa AbstractVector && length(output) == 4 &&
              all(value -> value isa Bool, output)
      valid || throw(ArgumentError(
        "Invalid categorifiability fixture $filename, case $index: " *
        "expected four Booleans [ZSC, Lagrange, PDC, D-number]. " *
        "Regenerate this fixture from the criterion results; do not reuse " *
        "the preceding universal-grading output variable.",
      ))
    end
  end

  return data
end

function load_anyonica_data(filename::AbstractString)
  path = anyonica_data_path(filename)
  isfile(path) || throw(ArgumentError("Missing Anyonica fixture: $path"))

  text = read(path, String)
  startswith(strip(text), "version https://git-lfs.github.com/spec/v1") &&
    throw(ArgumentError(
      "Fixture $path is a Git LFS pointer, not JSON data. " *
      "Retrieve the actual fixture with git lfs pull in the repository.",
    ))

  data = try
    JSON.parse(text)
  catch err
    err isa InterruptException && rethrow()
    throw(ArgumentError("Invalid JSON in fixture $path: " * sprint(showerror, err)))
  end

  return validate_anyonica_data(data, path)
end

# ============================================================
# JSON conversion and shape checks
# ============================================================

json_int(x::Bool) = throw(ArgumentError("A Boolean is not a fixture integer."))
json_int(x::Integer) = Int(x)

function json_int(x::AbstractFloat)
  isfinite(x) && isinteger(x) || throw(ArgumentError(
    "Expected an integer-valued JSON number; got $(repr(x)).",
  ))
  return Int(x)
end

json_int(x::AbstractString) = parse(Int, strip(x))

function _is_json_integer(x)
  x isa Bool && return false
  x isa Integer && return true
  x isa AbstractFloat && return isfinite(x) && isinteger(x)
  x isa AbstractString && return tryparse(Int, strip(x)) !== nothing
  return false
end

looks_like_int_vector(x) = x isa AbstractVector && all(_is_json_integer, x)

function looks_like_formal_code(x)
  looks_like_int_vector(x) && length(x) == 4 || return false
  r, m, nnsd, index = json_int.(x)
  return r >= 1 && m >= 1 && 0 <= nnsd < r && index >= 1
end

function looks_like_matrix(x)
  x isa AbstractVector && !isempty(x) || return false
  all(row -> row isa AbstractVector, x) || return false
  n = length(first(x))
  return n > 0 && all(row -> length(row) == n && looks_like_int_vector(row), x)
end

function looks_like_3tensor(x)
  x isa AbstractVector && !isempty(x) || return false
  r = length(x)
  return all(slice -> slice isa AbstractVector && length(slice) == r &&
    all(row -> looks_like_int_vector(row) && length(row) == r, slice), x)
end

function json_int_vector(x)
  looks_like_int_vector(x) || throw(ArgumentError("Expected an integer vector."))
  return Int[json_int(value) for value in x]
end

function json_int_matrix(x)
  looks_like_matrix(x) || throw(ArgumentError(
    "Expected a nonempty rectangular integer matrix in nested JSON form.",
  ))
  return [json_int(x[i][j]) for i in eachindex(x), j in eachindex(first(x))]
end

function json_int_3tensor(x)
  looks_like_3tensor(x) || throw(ArgumentError(
    "Expected a nonempty cubic integer tensor in nested JSON form.",
  ))
  r = length(x)
  return [json_int(x[i][j][k]) for i in 1:r, j in 1:r, k in 1:r]
end

# ============================================================
# Ring construction from Anyonica references
# ============================================================

function ring_from_anyonica_mt(x)
  return fusion_ring(json_int_3tensor(x); skip_check = true)
end

function code_from_json(x)
  looks_like_formal_code(x) || throw(ArgumentError(
    "Expected an AnyonWiki code [rank, multiplicity, nnsd, index]; got $(repr(x)).",
  ))
  return json_int_vector(x)
end

function ring_from_anyonica_code(x; fresh::Bool = false)
  stored = from_anyonwiki_code(code_from_json(x))
  fresh || return stored

  # Rebuild from the multiplication table, without cached invariants or UUID links.
  # This lets domain tests distinguish computation from catalogue retrieval.
  return fusion_ring(copy(multiplication_table(stored)); labels = copy(labels(stored)))
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
    error("Could not interpret Anyonica ring reference of type $(typeof(x)): $(repr(x))")
  end
end

# Alias with a short name for tests.
function ring_from_ref(x)
  return ring_from_anyonica_ref(x)
end

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

function check_mt_equal(actual_mt::Array{Int, 3}, expected_json, msg::AbstractString)
  expected = json_int_3tensor(expected_json)
  return check_equal_tensor(actual_mt, expected, msg)
end

function check_mt_equal(actual_ring::FusionRing, expected_json, msg::AbstractString)
  return check_mt_equal(multiplication_table(actual_ring), expected_json, msg)
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
  out = Tuple{Vector{Int}, Array{Int, 3}}[]

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
  out = Tuple{Vector{Int}, Array{Int, 3}}[]

  for item in xs
    if item isa Tuple
      els = collect(Int, item[1])
      subring = item[2]
      push!(out, (els, multiplication_table(subring)))
    elseif item isa Dict
      inj = if haskey(item, "injection")
        item["injection"]
      elseif haskey(item, :injection)
        item[:injection]
      else
        error("Subring Dict did not contain an injection key.")
      end

      sub = if haskey(item, "fusion_ring")
        item["fusion_ring"]
      elseif haskey(item, :fusion_ring)
        item[:fusion_ring]
      else
        error("Subring Dict did not contain a fusion_ring key.")
      end

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

# ============================================================
# Numeric reference comparisons
# ============================================================

function json_float(x::Union{Integer,AbstractFloat,AbstractString})
  x isa Bool && throw(ArgumentError("A Boolean is not a numeric reference value."))
  value = x isa AbstractString ? parse(Float64, strip(x)) : Float64(x)
  isfinite(value) || throw(ArgumentError("Expected a finite numeric reference value."))
  return value
end

json_float_vector(xs) = Float64[json_float(x) for x in xs]
json_nested_float_vector(xs) = [json_float_vector(row) for row in xs]
numeric_vector_from_json(xs) = json_float_vector(xs)

function check_numeric_vector_approx(
  actual, expected_json, msg::AbstractString; atol = 1e-8, rtol = 1e-8
)
  # Preserve imaginary parts: a large imaginary error must not disappear via real(x).
  actual_vec = ComplexF64[ComplexF64(x) for x in actual]
  expected_vec = numeric_vector_from_json(expected_json)
  return check_approx(actual_vec, expected_vec, msg; atol = atol, rtol = rtol)
end

function check_numeric_approx(
  actual, expected_json, msg::AbstractString; atol = 1e-8, rtol = 1e-8
)
  return check_approx(
    ComplexF64(actual), json_float(expected_json), msg; atol = atol, rtol = rtol
  )
end

# ============================================================
# Reference selection
# ============================================================

function oracle_mode()
  mode = lowercase(strip(get(ENV, "FUSIONRINGS_ORACLE_MODE", "representative")))
  mode in ("all", "full") && return :all
  mode == "first" && return :first
  mode == "spread" && return :spread
  mode in ("representative", "rep") && return :representative
  mode in ("small", "small_ranks") && return :small_ranks
  throw(ArgumentError("Unknown FUSIONRINGS_ORACLE_MODE=$(repr(mode))."))
end

function _rank_from_oracle_input(x)
  if looks_like_formal_code(x)
    return json_int(x[1])
  elseif looks_like_3tensor(x)
    return length(x)
  elseif x isa AbstractVector && length(x) == 2 && looks_like_formal_code(x[1])
    # Tensor-product inputs are [code1, code2]; their output rank is the product.
    if looks_like_formal_code(x[2])
      return json_int(x[1][1]) * json_int(x[2][1])
    end
    # Permutation inputs are [code, permutations]; rank is unchanged.
    return json_int(x[1][1])
  end
  return nothing
end

"""
    oracle_case_indices(data; mode=oracle_mode(), max_cases=12, max_rank=6)

first: at most max_cases initial cases.
spread / representative: at most max_cases evenly spaced cases.
all: every case.
small_ranks: every case whose input implies rank <= max_rank (not count-limited).

"representative" is a legacy name for positional sampling, not a guarantee of
mathematical coverage. Domain tests must add explicit cases for important branches.
"""
function oracle_case_indices(data; mode = oracle_mode(), max_cases = 12, max_rank = 6)
  mode in (:all, :first, :spread, :representative, :small_ranks) ||
    throw(ArgumentError("Unknown oracle selection mode: $mode"))
  max_cases isa Integer && !(max_cases isa Bool) && max_cases >= 1 ||
    throw(ArgumentError("max_cases must be a positive integer."))
  max_rank isa Integer && !(max_rank isa Bool) && max_rank >= 1 ||
    throw(ArgumentError("max_rank must be a positive integer."))

  n = length(data["Input"])
  n > 0 || throw(ArgumentError("Cannot select cases from an empty fixture."))

  if mode == :all
    return collect(1:n)
  elseif mode == :first
    return collect(1:min(max_cases, n))
  elseif mode in (:spread, :representative)
    count = min(max_cases, n)
    count == 1 && return [1]
    return sort(unique(round.(Int, range(1, n; length = count))))
  end

  indices = Int[]
  for index in 1:n
    r = _rank_from_oracle_input(data["Input"][index])
    r === nothing && throw(ArgumentError(
      "Cannot infer rank for oracle case $index. Use the filename-aware " *
      "selector for creation data, or choose a different mode.",
    ))
    r <= max_rank && push!(indices, index)
  end
  isempty(indices) && throw(ArgumentError("No oracle cases have rank <= $max_rank."))
  return indices
end

function oracle_case_indices(
  filename::AbstractString, data; mode = oracle_mode(), max_cases = 12, max_rank = 6
)
  creation_fixture = basename(filename) in (
    "zn_tables.json", "grouptables.json", "HI_tables.json", "TY_tables.json",
    "psu2k_tables.json", "su2k_tables.json", "son2_tables.json",
  )

  if creation_fixture && mode == :representative
    # Preserve the existing default: run all examples of these constructor families.
    return oracle_case_indices(data; mode = :all, max_cases = max_cases, max_rank = max_rank)
  elseif creation_fixture && mode == :small_ranks
    # The parameter is not always the rank (e.g. SU(2)_k has rank k+1).
    # Creation fixtures contain one raw output tensor per input.
    oracle_case_indices(data; mode = :first, max_cases = max_cases, max_rank = max_rank)
    indices = Int[]
    for index in eachindex(data["Output"])
      output = data["Output"][index]
      looks_like_3tensor(output) || throw(ArgumentError(
        "Invalid multiplication-table output in $filename, case $index.",
      ))
      length(output) <= max_rank && push!(indices, index)
    end
    isempty(indices) && throw(ArgumentError("No oracle cases have rank <= $max_rank in $filename."))
    return indices
  end

  return oracle_case_indices(data; mode = mode, max_cases = max_cases, max_rank = max_rank)
end
