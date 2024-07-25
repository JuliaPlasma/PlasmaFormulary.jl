using PlasmaFormulary
using Documenter

DocMeta.setdocmeta!(PlasmaFormulary, :DocTestSetup, :(using PlasmaFormulary); recursive=true)

makedocs(;
    modules=[PlasmaFormulary],
    authors="Luke Adams <luke@lukeclydeadams.com> and contributors",
    sitename="PlasmaFormulary.jl",
    format=Documenter.HTML(;
        canonical="https://adamslc.github.io/PlasmaFormulary.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/adamslc/PlasmaFormulary.jl",
    devbranch="main",
)
