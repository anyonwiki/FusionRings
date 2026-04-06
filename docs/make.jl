using FusionRings
using Documenter

DocMeta.setdocmeta!(FusionRings, :DocTestSetup, :(using FusionRings); recursive=true)

makedocs(;
    modules=[FusionRings],
    authors="gert-vercleyen <gert.vercleyen@protonmail.com>, Szagha02",
    sitename="FusionRings.jl",
    format=Documenter.HTML(;
        canonical="https://gert-vercleyen.gitlab.io/FusionRings.jl",
        edit_link="develop",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
