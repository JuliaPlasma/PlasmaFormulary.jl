```@meta
CurrentModule = PlasmaFormularyZ
```

# PlasmaFormularyZ

A Julia package for plasma physics formulas.
This package implements a subset of the formulas found in the [NRL Plasma Formulary](https://www.nrl.navy.mil/News-Media/Publications/NRL-Plasma-Formulary/), as well as some additional formulas taken from the [PlasmaPy](https://docs.plasmapy.org) project.

## Features

- For functions that take a `ParticleLike` as an argument, we use [ChargedParticles.jl](https://github.com/Beforerr/ChargedParticles.jl) package to handle particle properties. One could provide a symbolic particle name or string or a Particle object for the argument. In addition, mass number and charge number can be specified as keyword arguments for these functions.
- For function whose arguments are `Unitful` quantities with unique dimensions, arguments are order independent (see [PermuteArgs.jl](https://github.com/Beforerr/PermuteArgs.jl)).

## Quick Example

```@repl share
using PlasmaFormularyZ, Unitful

debye_length(1e18u"cm^-3", 10u"eV")

gyrofrequency(0.01u"T", :e) # electron gyrofrequency

plasma_frequency(1e19u"m^-3", "proton") # proton plasma frequency
```

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/Beforerr/PlasmaFormularyZ.jl")
```

## Related packages

- [PlasmaPy](https://docs.plasmapy.org) : `plasmapy.formulary` subpackage provides theoretical formulas for calculation of physical quantities helpful for plasma physics.
- [FusionFormulary.jl](https://github.com/JuliaFusion/FusionFormulary.jl)

## Index

```@index
```
