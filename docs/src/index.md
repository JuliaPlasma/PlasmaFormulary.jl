```@meta
CurrentModule = PlasmaFormulary
```

# PlasmaFormulary

Documentation for [PlasmaFormulary](https://github.com/Beforerr/PlasmaFormulary.jl).

We intensively use Julia's dynamic dispatch to handle arguments. As a result, arguments can be given in any order for most functions, since the type of each argument is often unique. See [@permutable_args](@ref) for more information.

For functions that take a `ParticleLike` as an argument, we use [ChargedParticles.jl](https://github.com/Beforerr/ChargedParticles.jl) to handle particle properties. One could provide a symbolic particle name or string or a `Particle` object for the argument. Also, mass and charge can be specified as keyword arguments for these functions. 

https://beforerr.github.io/ChargedParticles.jl/dev/

```@index
```