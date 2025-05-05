using PlasmaFormularyZ
using Documenter

DocMeta.setdocmeta!(
    PlasmaFormularyZ,
    :DocTestSetup,
    :(using PlasmaFormularyZ, Unitful, DimensionfulAngles);
    recursive = true,
)

makedocs(;
    modules = [PlasmaFormularyZ],
    sitename = "PlasmaFormularyZ.jl",
    pages = ["Home" => "index.md", "API Reference" => "api.md"],
    doctest = true,
)

deploydocs(;
    repo = "github.com/Beforerr/PlasmaFormularyZ.jl",
    devbranch = "main",
    push_preview = true,
)
