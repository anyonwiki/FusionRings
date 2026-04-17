using Test
using FusionRings

include("test_utils.jl")
include("test_fixtures.jl")

const TEST_GROUP = lowercase(get(ENV, "FUSIONRINGS_TEST_GROUP", "all"))
const TEST_LEVEL = lowercase(get(ENV, "FUSIONRINGS_TEST_LEVEL", "all"))

"""
    want_group(name)

Return true if the test group `name` should run.
"""
want_group(name::AbstractString) =
    TEST_GROUP == "all" || lowercase(name) == TEST_GROUP

"""
    want_level(name)

Return true if the test level `name` should run.
"""
want_level(name::AbstractString) =
    TEST_LEVEL == "all" || lowercase(name) == TEST_LEVEL

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
        include("subring_tests.jl")
    end

    if want_group("automorphism")
        include("automorphism_tests.jl")
    end

    if want_group("numeric")
        include("numeric_tests.jl")
    end

    if want_group("import_export")
        include("import_export_tests.jl")
    end

    if want_group("formatting")
        include("formatting_tests.jl")
    end
end