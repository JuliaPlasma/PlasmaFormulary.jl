using PlasmaFormulary
using Documenter

DocMeta.setdocmeta!(
    PlasmaFormulary,
    :DocTestSetup,
    :(using PlasmaFormulary, Unitful);
    recursive = true,
)

makedocs(;
    modules = [PlasmaFormulary],
    sitename = "PlasmaFormulary.jl",
    pages = ["Home" => "index.md", "API Reference" => "api.md"],
)

deploydocs(; repo = "github.com/JuliaPlasma/PlasmaFormulary.jl", devbranch = "main", push_preview = true)
