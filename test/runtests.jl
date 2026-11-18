using Test
using FusionRings

using FusionRings: fix_fractions, factor_squares, bold_integer

include("test_utils.jl")
include("test_fixtures.jl")

# ============================================================
# Configuration
# ============================================================

const TEST_GROUPS = (
  "creation",
  "properties",
  "operations",
  "subring",
  "automorphism",
  "numeric",
  "categorifiability",
  "import_export",
  "formatting",
)

const TEST_LEVELS = (
  "basic",
  "intermediate",
  "reference",
)

function validate_test_group(name::AbstractString)
  group = lowercase(strip(name))

  if group != "all" && !(group in TEST_GROUPS)
    throw(ArgumentError(
      "Unknown test group $(repr(name)). Expected one of: " *
      join(TEST_GROUPS, ", ") * ", all.",
    ))
  end

  return group
end

function validate_test_level(name::AbstractString)
  level = normalize_test_level(name)

  if level != "all" && !(level in TEST_LEVELS)
    throw(ArgumentError(
      "Unknown test level $(repr(name)). Expected basic, intermediate, " *
      "reference, all, or a supported alias such as level1.",
    ))
  end

  return level
end

const TEST_GROUP = validate_test_group(
  get(ENV, "FUSIONRINGS_TEST_GROUP", "all"),
)

const TEST_LEVEL = validate_test_level(
  get(ENV, "FUSIONRINGS_TEST_LEVEL", "all"),
)

# ============================================================
# Selection
# ============================================================

"""
    want_group(name)

Whether the requested group is enabled.
The group "all" imposes no additional restriction.
"""
function want_group(name::AbstractString)
  group = validate_test_group(name)

  return TEST_GROUP == "all" ||
         group == "all" ||
         group == TEST_GROUP
end

"""
    want_level(name)

Whether the requested level is enabled.
The level "all" imposes no additional restriction.
"""
function want_level(name::AbstractString)
  level = validate_test_level(name)

  return TEST_LEVEL == "all" ||
         level == "all" ||
         level == TEST_LEVEL
end

"""
    maybe_testset(name; group="all", level="all") do
        ...
    end

Run the body when both its group and level are enabled.
Prefer this explicit keyword form for new tests.
"""
function maybe_testset(
  f::Function,
  name::AbstractString;
  group::AbstractString = "all",
  level::AbstractString = "all",
)
  # Validate both arguments even when one selector excludes the test.
  group_enabled = want_group(group)
  level_enabled = want_level(level)

  if group_enabled && level_enabled
    @testset "$name" begin
      f()
    end
  end

  return nothing
end

"""
    maybe_testset(level, name) do
        ...
    end

Compatibility with the level-first calls used by the current test files.
"""
function maybe_testset(
  f::Function,
  level::AbstractString,
  name::AbstractString,
)
  return maybe_testset(f, name; level = level)
end

"""
    maybe_testset(name, group, level) do
        ...
    end

Compatibility with the previous three-string interface.
"""
function maybe_testset(
  f::Function,
  name::AbstractString,
  group::AbstractString,
  level::AbstractString,
)
  return maybe_testset(f, name; group = group, level = level)
end

# ============================================================
# File loading
# ============================================================

function include_test_file(filename::AbstractString)
  path = joinpath(@__DIR__, filename)

  isfile(path) ||
    error("Selected test file is missing: $path")

  # The current subring, automorphism, and numeric files are blank.
  # Record these as skipped coverage rather than silently doing nothing.
  if isempty(strip(read(path, String)))
    @info "Skipping empty test file" file = filename

    @testset "Empty test file: $filename" begin
      @test_skip false
    end

    return nothing
  end

  include(path)
  return nothing
end

# ============================================================
# Run selected groups
# ============================================================

@info "FusionRings test selection" group = TEST_GROUP level = TEST_LEVEL

@testset "FusionRings.jl" begin
  for group in TEST_GROUPS
    if want_group(group)
      include_test_file(group * "_tests.jl")
    end
  end
end