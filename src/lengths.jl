"""
Calculate the radius of circular motion for a charged particle in a uniform magnetic field

References: [PlasmaPy API Documentation](https://docs.plasmapy.org/en/latest/api/plasmapy.formulary.lengths.gyroradius.html)

# Examples
```jldoctest
julia> gyroradius(0.2u"T", 1e6u"K")
0.006682528174870854 m
```
"""
function gyroradius(B, mass, q, Vperp::Velocity)
    return upreferred(mass * Vperp / (q * B))
end

function gyroradius(B, Vperp::Velocity; mass_numb = 1, Z = 1)
    return gyroradius(B, mass_numb * Unitful.u, Z * Unitful.q, Vperp)
end

function gyroradius(B, T::EnergyOrTemp; mass_numb = 1, Z = 1)
    Vperp = thermal_speed(T, mass_numb)
    return gyroradius(B, Vperp; mass_numb, Z)
end

gyroradius(val::Union{Velocity,EnergyOrTemp}, B::BField; kw...) = gyroradius(B, val; kw...)

electron_gyroradius(B, Vperp::Velocity) = gyroradius(B, me, Unitful.q, Vperp)
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
The inertial length is the characteristic length scale for a particle to be accelerated in a plasma. The Hall effect becomes important on length scales shorter than the ion inertial length.

References: [PlasmaPy API Documentation](https://docs.plasmapy.org/en/latest/api/plasmapy.formulary.lengths.inertial_length.html)
"""
function inertial_length(n::NumberDensity, q, mass::Mass)
    upreferred(c / plasma_frequency(n, q, mass))
end

electron_inertial_length(n::NumberDensity) = upreferred(c / plasma_frequency(n))

function ion_inertial_length(n::NumberDensity, Z, mass::Mass)
    inertial_length(n, Z * Unitful.q, mass)
end

ion_inertial_length(n::NumberDensity; Z = 1, mass = mp) = ion_inertial_length(n, Z, mass)

function debye_length(density::NumberDensity, eot::EnergyOrTemp)
    upreferred(sqrt(ε0 * energy(eot) / density / q^2))
end