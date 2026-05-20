using Test
using FusionRings

include("test_utils.jl")
include("test_fixtures.jl")

# ============================================================
# Test group / level configuration
# ============================================================

const TEST_GROUP = lowercase(get(ENV, "FUSIONRINGS_TEST_GROUP", "all"))

# test_utils.jl defines normalize_test_level(...)
const TEST_LEVEL = normalize_test_level(get(ENV, "FUSIONRINGS_TEST_LEVEL", "all"))

"""
    want_group(name)

Return true if the test group `name` should run.

Examples:
    FUSIONRINGS_TEST_GROUP=creation
    FUSIONRINGS_TEST_GROUP=properties
    FUSIONRINGS_TEST_GROUP=all
"""
want_group(name::AbstractString) =
    TEST_GROUP == "all" || lowercase(name) == TEST_GROUP

"""
    want_level(name)

Return true if the test level `name` should run.

Accepted names:
    basic, intermediate, reference
    level1, level2, level3
    all
"""
want_level(name::AbstractString) =
    TEST_LEVEL == "all" || normalize_test_level(name) == TEST_LEVEL

"""
    include_if_exists(path)

Include  test file only if  exists.  This  useful for optional
test groups which may be empty or absent while the main tests
live in properties_tests.jl.
"""
function include_if_exists(path::AbstractString)
    fullpath = joinpath(@__DIR__, path)

    if isfile(fullpath)
        include(path)
    else
        @info "Skipping missing optional test file" path
    end
end

# ============================================================
# Main tests
# ============================================================

@testset "FusionRings.jl" begin

    if want_group("creation")
        include("creation_tests.jl")
    end

    if want_group("properties")
        include("properties_tests.jl")
    end

    if want_group("operations")
        include("operations_tests.jl")
    end

    if want_group("subring")
        include_if_exists("subring_tests.jl")
    end

    if want_group("automorphism")
        include_if_exists("automorphism_tests.jl")
    end

    if want_group("numeric")
        include_if_exists("numeric_tests.jl")
    end

    if want_group("import_export")
        include("import_export_tests.jl")
    end

    if want_group("formatting")
        include("formatting_tests.jl")
    end
end

#example: julia --project=. test/runtests.jl
#FUSIONRINGS_TEST_LEVEL=level1 julia --project=. test/runtests.jl
#FUSIONRINGS_TEST_LEVEL=level3 julia --project=. test/runtests.jl - representative anyonica tests
#FUSIONRINGS_TEST_LEVEL=level3 FUSIONRINGS_ORACLE_MODE=all julia --project=. test/runtests.jl -full anyonica tests
#only creation tests: FUSIONRINGS_TEST_GROUP=creation FUSIONRINGS_TEST_LEVEL=level3 julia --project=. test/runtests.jl
#level1 = basic hand-written sanity tests
#level2 = intermediate hand-checkable correctness tests
#level3 = reference / Anyonica oracle tests