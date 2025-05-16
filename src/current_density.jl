"""
    Alfven_current_density(B::BFieldOrBFields, n::NumberDensity)

Calculate the Alfvén current density ``J_A``, a natural scaling for current density: 

```math
J_A = e n V_A = \\frac{B}{μ₀ dᵢ}
```

where ``V_A`` is the Alfven speed, ``dᵢ`` is ion inertial length.

See also: [`Alfven_speed`](@ref), [`inertial_length`](@ref).
"""
function Alfven_current_density end

@permute_args Alfven_current_density(Va::Velocity, n::NumberDensity) =
    q * n * Va |> upreferred
@permute_args Alfven_current_density(𝐁::BFieldOrBFields, di::Length) =
    norm(𝐁) / (μ0 * di) |> upreferred
@permute_args Alfven_current_density(𝐁::BFieldOrBFields, n::NumberDensity) =
    Alfven_current_density(Alfven_speed(𝐁, n), n)
