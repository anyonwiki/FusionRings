# Script to update Artifacts.toml file
# !!! Need to have git installed and web access for this to work !!!

using Pkg.Artifacts
using Pkg.PlatformEngines: unpack
using Downloads
using SHA

# path to artifacts toml in project root
atfn = joinpath( @__DIR__ , "Artifacts.toml" )


# when running script we want the data to come from newest commit
repo   = "anyonwiki/AnyonWikiDatabase"
branch = "main"
ref = ( String ∘ first ∘ split )( 
    read(
        `git ls-remote https://github.com/$repo.git refs/heads/$branch`, 
        String
    )
)

tarballs = [
    "FusionRings"      => "AlgebraicStructures/FusionRings.tar.gz",
    "AlgebraicNumbers" => "SupportingData/AlgebraicNumbers.tar.gz"
]

function bind_tarball!( atfn, name, pth )
    url = "https://raw.githubusercontent.com/$repo/$ref/$pth"

    mktempdir() do tmp 
        tarball = joinpath( tmp, basename( pth) )
        @info "Downloading $name" url
        Downloads.download(url,tarball)    

        tb_hash = (bytes2hex ∘ sha256 ∘ open)(tarball)
        
        up(dir)   = unpack( tarball, dir ) 
        tree_hash = create_artifact(up)

        bind_artifact!(
            atfn, name, tree_hash;
            download_info = [ ( url, tb_hash ) ],
            force = true
        )
        @info "Bound $name" tree_hash tb_hash
    end
end

for (name,path) in tarballs
    bind_tarball!(atfn,name,path)
end