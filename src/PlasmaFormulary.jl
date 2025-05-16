module PlasmaFormulary

using Unitful
using Unitful: 𝐋, 𝐈, 𝐌, 𝐓
using Unitful: μ0, ε0, c, q
using Unitful: k, ħ
using Unitful: me, mp, u
using Unitful: Velocity, Mass, Force, BField, EField, Density, Charge, Length
using UnitfulEquivalences
using DimensionfulAngles: radᵃ as rad
using PermuteArgs
using ChargedParticles
using ChargedParticles: ParticleLike
using LinearAlgebra: norm, norm2, ×, ⋅

@derived_dimension NumberDensity 𝐋^-3
@derived_dimension CurrentDensity 𝐈 / 𝐋^2
@derived_dimension PressureGradient 𝐌 * 𝐓^-2 / 𝐋^2
const EnergyOrTemp = Union{Unitful.Temperature,Unitful.Energy}
const BFields = AbstractVector{<:BField}
const BFieldOrBFields = Union{BField,BFields}
const EFields = AbstractVector{<:EField}
const Forces = AbstractVector{<:Force}
const PressureGradients = AbstractVector{<:PressureGradient}
const qe = Unitful.q

# Workaround for UnitfulEquivalences.uconvert not supporting conversions of
# quantites with the same dimensions. See
# https://github.com/sostock/UnitfulEquivalences.jl/issues/19
function custom_uconvert(dest_unit, x, equivalence)
    if dimension(x) == dimension(dest_unit)
        # plain Unitful ustrip
        Unitful.uconvert(dest_unit, x)
    else
        UnitfulEquivalences.uconvert(dest_unit, x, equivalence)
    end
end

energy(eot) = custom_uconvert(u"J", eot, Thermal())
temperature(eot) = custom_uconvert(u"K", eot, Thermal())

include("constants.jl")
include("dimensionless.jl")
export plasma_beta

include("lengths.jl")
const debye_length = Debye_length # to be removed
export gyroradius,
    electron_gyroradius,
    electron_debroglie_length,
    classical_minimum_approach_distance,
    inertial_length,
    electron_inertial_length,
    ion_inertial_length,
    Debye_length,
    debye_length

include("speeds.jl")
include("velocities.jl")
export Alfven_velocity, Alfven_speed, ion_sound_speed, thermal_speed
export diamagnetic_drift, ExB_drift, force_drift
export thermal_temperature

include("current_density.jl")
export Alfven_current_density

include("frequencies.jl")
export gyrofrequency, plasma_frequency

include("misc.jl")
export thermal_pressure, magnetic_pressure

include("alias.jl")

end
