module PlasmaFormulary

using Unitful
using UnitfulEquivalences
using Unitful: μ0, ε0, c, q
using Unitful: k, ħ
using Unitful: me, mp, u
using Unitful: Velocity, Mass, BField, Density, Charge
using LinearAlgebra
using Combinatorics
using ChargedParticles: Particle, ParticleLike

@derived_dimension NumberDensity Unitful.𝐋^-3

const EnergyOrTemp = Union{Unitful.Temperature,Unitful.Energy}
energy(eot) = uconvert(u"J", eot, Thermal())
temperature(eot) = uconvert(u"K", eot, Thermal())
temperature(T::Unitful.Temperature) = uconvert(u"K", T)

# Physical Constants (SI)
rydberg_constant = 1.09737e7 * u"m^-1"
bohr_radius = 5.29177e-11 * u"m"
atomic_cross_section = 8.79736e-21 * u"m^2"
classical_electron_radius = 2.81794e-15 * u"m"
thomson_cross_section = 6.65246e-29 * u"m^2"
electron_compton_wavelength = 2.42631e-12 * u"m"
reduced_electron_compton_wavelength = 3.86159e-13 * u"m"
fine_structure_constant = 7.29735e-3
# Scales associated with 1 eV excluded
# Energy scales excluded
avogadro_number = 6.02214e23 * u"mol^-1"
faraday_constant = 9.64853e4 * u"C / mol"
gas_constant = 8.31446 * u"J / K / mol"
loschmidt_number = 2.68678e25 * u"m^-3"
atomic_mass_unit = 1.66054e-27 * u"kg"
standard_temperature = 273.15 * u"K"
atmospheric_pressure = 1.01325e5 * u"Pa"
# Pressure of 1 torr excluded
standard_molar_volume = 2.24140e-2 * u"m^3 / mol"
molar_weight_of_air = 2.89647e-2 * u"kg / mol"

# Fundamental plasma parameters
# These formulas have been converted to use SI units from the original Gaussian cgs units
# that are used in the 2023 edition of the formulary.
include("utils.jl")
include("dimensionless.jl")
include("lengths.jl")
include("speeds.jl")
include("frequencies.jl")
include("misc.jl")

export plasma_beta
export gyroradius, debye_length, inertial_length
export Alfven_velocity, Alfven_speed
export gyrofrequency, plasma_frequency

# aliases
const rc_ = gyroradius
export rc_

end
