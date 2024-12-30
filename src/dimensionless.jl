"""
    plasma_beta(T, n, B)

Compute the plamsa beta (β), the ratio of thermal pressure to magnetic pressure.

# Arguments
- `T`: The temperature of the plasma.
- `n`: The particle density of the plasma.
- `B`: The magnetic field in the plasma.

# Example
```jldoctest
julia> plasma_beta(1e6u"K", 1e19u"m^-3", 0.2u"T")
0.008674873511172188
```

# References
- [PlasmaPy Documentation](https://docs.plasmapy.org/en/stable/api/plasmapy.formulary.dimensionless.beta.html)
- [Wikipedia](https://en.wikipedia.org/wiki/Plasma_beta)
"""
function plasma_beta end
@permute_args function plasma_beta(T::EnergyOrTemp, n::NumberDensity, B::Unitful.BField)
    p_th = thermal_pressure(T, n)
    p_B = magnetic_pressure(B)
    return p_th / p_B |> upreferred
end
