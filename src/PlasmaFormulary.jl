module PlasmaFormulary

using Unitful
using UnitfulEquivalences
using Unitful: μ0, ε0, c, q
using Unitful: k, ħ
using Unitful: me, mp, u
using Unitful: Velocity, Mass, BField, Density, Charge

@derived_dimension NumberDensity Unitful.𝐋^-3

const EnergyOrTemp = Union{Unitful.Temperature, Unitful.Energy}

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

include("dimensionless.jl")
export plasma_beta

include("lengths.jl")
export gyroradius, electron_gyroradius, electron_debroglie_length, classical_minimum_approach_distance, inertial_length, electron_inertial_length, ion_inertial_length, debye_length

include("speeds.jl")
export alfven_velocity, thermal_velocity, thermal_temperature, electron_thermal_velocity

include("frequencies.jl")
export gyrofrequency, electron_gyrofrequency, ion_gyrofrequency, plasma_frequency

include("misc.jl")
export thermal_pressure, magnetic_pressure


end
