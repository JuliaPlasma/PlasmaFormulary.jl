"""
    gyrofrequency(B::BField, p::ParticleLike; kw...)
    gyrofrequency(B::BField, mass::Mass, q::Charge; to_hz = false)

Calculate the gyrofrequency (or cyclotron frequency) of a charged particle's circular motion in a magnetic field.
The gyrofrequency is the angular frequency of a charged particle's gyromotion around magnetic field lines.

Set `to_hz = true` to convert output from angular frequency to Hz.

References: 
- [PlasmaPy Documentation](https://docs.plasmapy.org/en/latest/api/plasmapy.formulary.frequencies.gyrofrequency.html)

# Examples
```jldoctest; filter = r"(\\^-1|⁻¹)"
julia> gyrofrequency(0.01u"T", :p)  # proton gyrofrequency
957883.3292211705 rad s⁻¹

julia> gyrofrequency(0.1u"T", :e; to_hz = true)  # electron gyrofrequency as frequency
2.799248987233304e9 Hz

julia> gyrofrequency(250u"Gauss", "Fe"; z=13)  # Fe2+ ion gyrofrequency
560682.3520611608 rad s⁻¹
```
"""
function gyrofrequency end

@permute_args function gyrofrequency(B::BField, mass::Mass, q::Charge; to_hz = false)
    ω = upreferred(rad * abs(q * B / mass))
    return _to_Hz(ω, to_hz)
end

@permute_args function gyrofrequency(B::BField, p::ParticleLike; to_hz = false, kw...)
    p = particle(p; kw...)
    return gyrofrequency(B, mass(p), charge(p); to_hz)
end

"""
    plasma_frequency(n::NumberDensity, [q::Charge, mass::Mass]; to_hz = false)
    plasma_frequency(n::NumberDensity, p::ParticleLike; kw...)

Calculate the plasma frequency of a species.

```math
ω_p = \\sqrt{\\frac{n q^2}{m ε₀}}
```

The plasma frequency is a characteristic frequency of the plasma. 
More often, it refers to the frequency at which electrons oscillate in the plasma.

Set `to_hz = true` to convert output from angular frequency to Hz.

# Examples
```jldoctest; filter = r"(\\^-1|⁻¹)"
julia> plasma_frequency(1e19u"m^-3")  # plasma frequency
1.7839863654934653e11 rad s⁻¹

julia> plasma_frequency(1e19u"m^-3", :p)  # proton plasma frequency
4.1632945624883513e9 rad s⁻¹
```

# References
- [Wikipedia](https://en.wikipedia.org/wiki/Plasma_oscillation)
- [PlasmaPy Documentation](https://docs.plasmapy.org/en/latest/api/plasmapy.formulary.frequencies.plasma_frequency.html)
"""
function plasma_frequency end

@permute_args function plasma_frequency(n::NumberDensity, q::Charge, mass::Mass; to_hz = false)
    ω = rad * abs(q) * sqrt(n / mass / ε0) |> upreferred
    return _to_Hz(ω, to_hz)
end

@permute_args function plasma_frequency(n::NumberDensity, p::ParticleLike; to_hz = false, kw...)
    p = particle(p; kw...)
    return plasma_frequency(n, charge(p), mass(p); to_hz)
end

plasma_frequency(n::NumberDensity; to_hz = false) = plasma_frequency(n, Unitful.q, me; to_hz)

"""
    lower_hybrid_frequency(B::BField, n_i::NumberDensity, ion::ParticleLike; kw...)

Calculate the lower hybrid frequency.

```math
\\frac{1}{ω_{lh}^2} = \\frac{1}{ω_{ci}^2 + ω_{pi}^2} + \\frac{1}{ω_{ci} ω_{ce}}
```

where ω_ci is the ion gyrofrequency, ω_ce is the electron gyrofrequency,
and ω_pi is the ion plasma frequency.

Set `to_hz = true` to convert output from angular frequency to Hz.

# Examples
```jldoctest; filter = r"(\\^-1|⁻¹)"
julia> lower_hybrid_frequency(0.2u"T", 5e19u"m^-3", "D+")
5.783727350182709e8 rad s⁻¹
```

# References
- [Wikipedia](https://en.wikipedia.org/wiki/Lower_hybrid_oscillation)
- [PlasmaPy Documentation](https://docs.plasmapy.org/en/latest/api/plasmapy.formulary.frequencies.lower_hybrid_frequency.html)
"""
function lower_hybrid_frequency end

@permute_args function lower_hybrid_frequency(
    B::BField,
    n_i::NumberDensity,
    ion::ParticleLike;
    to_hz = false,
    kw...,
)
    p = particle(ion; kw...)
    ω_ci = gyrofrequency(B, mass(p), charge(p))
    ω_ce = gyrofrequency(B, me, Unitful.q)
    ω_pi = plasma_frequency(n_i, charge(p), mass(p))

    ω_lh = sqrt(inv(inv(ω_ci^2 + ω_pi^2) + inv(ω_ci * ω_ce)))
    return _to_Hz(ω_lh, to_hz)
end

"""
    upper_hybrid_frequency(B::BField, n_e::NumberDensity; to_hz = false)

Calculate the upper hybrid frequency.

```math
ω_{uh}^2 = ω_{ce}^2 + ω_{pe}^2
```

where ω_ce is the electron gyrofrequency and ω_pe is the electron plasma frequency.

Set `to_hz = true` to convert output from angular frequency to Hz.

# Examples
```jldoctest; filter = r"(\\^-1|⁻¹)"
julia> upper_hybrid_frequency(0.2u"T", 5e19u"m^-3")
4.004594195988481e11 rad s⁻¹

julia> upper_hybrid_frequency(0.2u"T", 5e19u"m^-3"; to_hz = true)
6.37350961368681e10 Hz
```

# References
- [Wikipedia](https://en.wikipedia.org/wiki/Upper_hybrid_oscillation)
- [PlasmaPy Documentation](https://docs.plasmapy.org/en/latest/api/plasmapy.formulary.frequencies.upper_hybrid_frequency.html)
"""
function upper_hybrid_frequency end

@permute_args function upper_hybrid_frequency(B::BField, n_e::NumberDensity; to_hz = false)
    ω_ce = gyrofrequency(B, me, Unitful.q)
    ω_pe = plasma_frequency(n_e)

    ω_uh = sqrt(ω_ce^2 + ω_pe^2)
    return _to_Hz(ω_uh, to_hz)
end
