using PlasmaFormulary
using Documenter

DocMeta.setdocmeta!(
    PlasmaFormulary,
    :DocTestSetup,
    :(using PlasmaFormulary; using Unitful);
    recursive = true,
)

makedocs(;
    modules = [PlasmaFormulary],
    sitename = "PlasmaFormulary.jl",
    pages = ["Home" => "index.md", "API Reference" => "api.md"],
)

deploydocs(; repo = "github.com/Beforerr/PlasmaFormulary.jl", devbranch = "main")
