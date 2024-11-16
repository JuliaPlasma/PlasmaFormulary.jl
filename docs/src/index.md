```@meta
CurrentModule = PlasmaFormulary
```

# PlasmaFormulary

A Julia package for plasma physics formulas. This package implements a subset of the formulas found in the [NRL Plasma Formulary](https://www.nrl.navy.mil/News-Media/Publications/NRL-Plasma-Formulary/).

We use Julia's dynamic dispatch to handle arguments. As a result, arguments can be given in any order for most functions, since the unit (type) of each argument is often unique. See [@permutable_args](@ref) for implementations of permutable arguments.

For functions that take a `ParticleLike` as an argument, we use [ChargedParticles.jl](https://github.com/Beforerr/ChargedParticles.jl) to handle particle properties. One could provide a symbolic particle name or string or a `Particle` object for the argument. Also, mass number and charge number can be specified as keyword arguments for these functions.

## Quick Example

```@repl share
using PlasmaFormulary, Unitful

gyrofrequency(0.1u"T", :e)  # electron gyrofrequency

plasma_frequency(1e19u"m^-3")  # electron plasma frequency

gyroradius(0.1u"T", :p, 1e6u"K")  # proton gyroradius at 1 million K
```

One does not need to care the order of the arguments, we automatically handle them for you. 🥳

Common particle names are also supported. 🫡

```@repl share
gyroradius(0.1u"T", :e, 1e6u"K") == 
    gyroradius(1e6u"K", "electron", 0.1u"T") == 
    gyroradius("e+", 1e6u"K", 0.1u"T")
```

## Related packages

- [PlasmaPy](https://docs.plasmapy.org) : `plasmapy.formulary` subpackage provides theoretical formulas for calculation of physical quantities helpful for plasma physics.

## Index

```@index
```