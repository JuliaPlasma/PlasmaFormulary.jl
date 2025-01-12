```@meta
CurrentModule = PlasmaFormulary
```

# PlasmaFormulary

A Julia package for plasma physics formulas.
This package implements a subset of the formulas found in the [NRL Plasma Formulary](https://www.nrl.navy.mil/News-Media/Publications/NRL-Plasma-Formulary/), as well as some additional formulas taken from the [PlasmaPy](https://docs.plasmapy.org) project.

In the future, the package will support providing particle properties using the [ChargedParticles.jl](https://github.com/Beforerr/ChargedParticles.jl) package.

## Quick Example

```@repl share
using PlasmaFormulary, Unitful

PlasmaFormulary.debye_length(1e18u"cm^-3", 10u"eV")
```

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/JuliaPlasma/PlasmaFormulary.jl")
```

## Related packages

- [PlasmaPy](https://docs.plasmapy.org) : `plasmapy.formulary` subpackage provides theoretical formulas for calculation of physical quantities helpful for plasma physics.
- [FusionFormulary.jl](https://github.com/JuliaFusion/FusionFormulary.jl)

## Index

```@index
```
