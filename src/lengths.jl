"""
    gyroradius(B::BField, p::ParticleLike, Vperp::Velocity; kw...)
    gyroradius(B::BField, p::ParticleLike, T::EnergyOrTemp; kw...)

Calculate the radius of circular motion for a charged particle in a uniform magnetic field.

Internal function:

    _gyroradius(B::BField, mass::Mass, q::Charge, Vperp::Velocity)

References: [PlasmaPy API Documentation](https://docs.plasmapy.org/en/latest/api/plasmapy.formulary.lengths.gyroradius.html)

# Examples
```jldoctest
julia> gyroradius(0.2u"T", :p, 1e6u"K")  # proton at 1 million K
0.006706796656500042 m

julia> gyroradius(0.2u"T", :e, 1e6u"K")  # electron at 1 million K
0.00015651672339994665 m

julia> gyroradius(0.2u"T", "Fe2+", 1e6u"K")
0.024988956222685512 m
```
"""
function gyroradius end

_gyroradius(B::BField, mass::Mass, q::Charge, Vperp::Velocity) =
    upreferred(abs(mass * Vperp / (q * B)))

@permutable_args function gyroradius(B::BField, p::ParticleLike, Vperp::Velocity; kw...)
    p = Particle(p; kw...)
    return _gyroradius(B, p.mass, p.charge, Vperp)
end

@permutable_args function gyroradius(B::BField, p::ParticleLike, T::EnergyOrTemp; kw...)
    p = Particle(p; kw...)
    Vperp = thermal_speed(T, p.mass)
    return _gyroradius(B, p.mass, p.charge, Vperp)
end

electron_gyroradius(B, Vperp::Velocity) = _gyroradius(B, me, Unitful.q, Vperp)
electron_gyroradius(B, T::EnergyOrTemp) = electron_gyroradius(B, thermal_speed(T, me))

# electron and ion trapping rates excluded

#electron and ion collision rates excluded

function electron_debroglie_length(eot::EnergyOrTemp)
    upreferred(sqrt(2 * pi * ħ^2 / me / energy(eot)))
end

function classical_minimum_approach_distance(eot::EnergyOrTemp)
    upreferred(q^2 / energy(eot) / (4 * pi * ε0))
end

"""
    inertial_length(n::NumberDensity, p::ParticleLike; kw...)

Calculate a charged particle's inertial length. 

The inertial length is the characteristic length scale for a particle to be accelerated in a plasma.
The Hall effect becomes important on length scales shorter than the ion inertial length.

References: [PlasmaPy API Documentation](https://docs.plasmapy.org/en/latest/api/plasmapy.formulary.lengths.inertial_length.html)
"""
function inertial_length end

_inertial_length(n::NumberDensity, q::Charge, mass::Mass) =
    upreferred(c / _plasma_frequency(n, q, mass))

@permutable_args function inertial_length(n::NumberDensity, p::ParticleLike; kw...)
    p = Particle(p; kw...)
    return _inertial_length(n, p.charge, p.mass)
end

electron_inertial_length(n::NumberDensity) = upreferred(c / plasma_frequency(n))


function debye_length(density::NumberDensity, eot::EnergyOrTemp)
    upreferred(sqrt(ε0 * energy(eot) / density / q^2))
end