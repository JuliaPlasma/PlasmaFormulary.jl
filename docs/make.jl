using PlasmaFormulary
using Documenter

DocMeta.setdocmeta!(
    PlasmaFormulary,
    :DocTestSetup,
    :(using PlasmaFormulary, Unitful, DimensionfulAngles);
    recursive = true,
)

makedocs(;
    modules = [PlasmaFormulary],
    sitename = "PlasmaFormulary.jl",
    pages = ["Home" => "index.md", "API Reference" => "api.md"],
    doctest = true,
)

deploydocs(;
    repo = "github.com/JuliaPlasma/PlasmaFormulary.jl",
    devbranch = "main",
    push_preview = true,
)
