using Pkg.Artifacts

# Path to the Artifacts.toml file in your package root
artifacts_toml = joinpath(dirname(@__DIR__), "Artifacts.toml")

# The folder containing the FusionRings data
data_source = joinpath(@__DIR__, "src", "data", "FusionRings" )

# Create the artifact from the data_source folder
FusionRings_hash = create_artifact() do artifact_dir
    # Copy or create your data files inside the artifact directory
    cp(data_source, artifact_dir; force=true)
end

# Define where the artifact tarball will be hosted (e.g., on GitHub Releases)
tarball_url = "https://github.com/anyonwiki/AnyonWikiDatabase/blob/main/FusionRings.tar.gz"
tarball_hash = archive_artifact(FusionRings_hash, "FusionRings.tar.gz")

# Bind the artifact to a name in Artifacts.toml
bind_artifact!(artifacts_toml, "FusionRings", FusionRings_hash;
               download_info = [(tarball_url, tarball_hash)],
               force = true)

# The folder containing the data on algebraic numbers
data_source = joinpath(@__DIR__, "src", "data", "AlgebraicNumbers" )

# Create the artifact from the data_source folder
AlgebraicNumbers_hash = create_artifact() do artifact_dir
    # Copy or create your data files inside the artifact directory
    cp(data_source, artifact_dir; force=true)
end

# Define where the artifact tarball will be hosted (e.g., on GitHub Releases)
tarball_url = "https://github.com/anyonwiki/AnyonWikiDatabase/blob/main/AlgebraicNumbers.tar.gz"
tarball_hash = archive_artifact(AlgebraicNumbers_hash, "AlgebraicNumbers.tar.gz")

# Bind the artifact to a name in Artifacts.toml
bind_artifact!(artifacts_toml, "AlgebraicNumbers", AlgebraicNumbers_hash;
               download_info = [(tarball_url, tarball_hash)],
               force = true)
