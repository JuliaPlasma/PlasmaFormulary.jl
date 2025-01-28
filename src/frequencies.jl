"""
    gyrofrequency(B::BField, p::ParticleLike; kw...)
    gyrofrequency(B::BField, mass::Mass, q::Charge)

Calculate the gyrofrequency (or cyclotron frequency) of a charged particle's circular motion in a magnetic field.
The gyrofrequency is the frequency of a charged particle's gyromotion around magnetic field lines.

References: 
- [PlasmaPy Documentation](https://docs.plasmapy.org/en/latest/api/plasmapy.formulary.frequencies.gyrofrequency.html)

# Examples
```jldoctest; filter = r"(\\^-1|⁻¹)"
julia> gyrofrequency(0.01u"T", :p)  # proton gyrofrequency
957883.3292211705 s⁻¹

julia> gyrofrequency(0.1u"T", :e)  # electron gyrofrequency
1.7588200107721634e10 s⁻¹

julia> gyrofrequency(250u"Gauss", "Fe"; z=13)  # Fe2+ ion gyrofrequency
560682.3520611608 s⁻¹
```
"""
function gyrofrequency end

@permute_args function gyrofrequency(B::BField, mass::Mass, q::Charge)
    upreferred(abs(q * B / mass))
end

@permute_args function gyrofrequency(B::BField, p::ParticleLike; kw...)
    p = particle(p; kw...)
    return gyrofrequency(B, mass(p), charge(p))
end

"""
    plasma_frequency(n::NumberDensity, [q::Charge, mass::Mass])
    plasma_frequency(n::NumberDensity, p::ParticleLike; kw...)

Calculate the plasma frequency of a species.

The plasma frequency is a characteristic frequency of the plasma. 
More often, it refers to the frequency at which electrons oscillate in the plasma.

# Examples
```jldoctest; filter = r"(\\^-1|⁻¹)"
julia> plasma_frequency(1e19u"m^-3")  # plasma frequency
1.7839863654934653e11 s⁻¹

julia> plasma_frequency(1e19u"m^-3", :p)  # proton plasma frequency
4.1632945624883513e9 s⁻¹
```

# References
- [Wikipedia](https://en.wikipedia.org/wiki/Plasma_oscillation)
- [PlasmaPy Documentation](https://docs.plasmapy.org/en/latest/api/plasmapy.formulary.frequencies.plasma_frequency.html)
"""
function plasma_frequency end

@permute_args function plasma_frequency(n::NumberDensity, q::Charge, mass::Mass)
    abs(q) * sqrt(n / mass / ε0) |> upreferred
end

@permute_args function plasma_frequency(n::NumberDensity, p::ParticleLike; kw...)
    p = particle(p; kw...)
    return plasma_frequency(n, charge(p), mass(p))
end

plasma_frequency(n::NumberDensity) = plasma_frequency(n, Unitful.q, me)
