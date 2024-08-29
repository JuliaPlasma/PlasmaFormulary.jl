module PlasmaFormulary

using Unitful

@derived_dimension NumberDensity Unitful.𝐋^-3

# Physical Constants (SI)
boltzmann_constant = 1.38065e-23 * u"J / K"
elementry_charge = 1.60219e-19 * u"C"
electron_mass = 9.10938e-31 * u"kg"
proton_mass = 1.67262e-27 * u"kg"
gravitational_constant = 6.67430e-11 * u"m^3 / s^2 / kg"
planck_constant = 6.62607e-34 * u"J * s"
reduced_planck_constant = 1.05457e-34 * u"J * s"
light_speed = 2.99792e8 * u"m / s"
free_space_permittivity = 8.85419e-12 * u"F / m"
free_space_permeability = 1.25662e-6 *u"H / m"
# Proton/electron mass ratio excluded
# Electron charge/mass ratio excluded
rydberg_constant = 1.09737e7 * u"m^-1"
bohr_radius = 5.29177e-11 * u"m"
atomic_cross_section = 8.79736e-21 * u"m^2"
classical_electron_radius = 2.81794e-15 * u"m"
thomson_cross_section = 6.65246e-29 * u"m^2"
electron_compton_wavelength = 2.42631e-12 * u"m"
reduced_electron_compton_wavelength = 3.86159e-13 * u"m"
fine_structure_constant = 7.29735e-3
# Radiation constants excluded
stefan_boltzmann_constant = 5.67037e-8 * u"W / m^2 / K^4"
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
# Calorie conversion omitted
gravitational_acceleration = 9.80665 * u"m / s^2"

# Fundamental plasma parameters
# These formulas have been converted to use SI units from the original Gaussian cgs units
# that are used in the 2023 edition of the formulary.
function electron_gyrofrequency(magnetic_field::Unitful.BField)
    upreferred(elementry_charge * magnetic_field / electron_mass)
end

function ion_gyrofrequency(magnetic_field::Unitful.BField, charge_state::Integer, ion_mass::Unitful.Mass)
    upreferred(charge_state * elementry_charge * magnetic_field / ion_mass)
end

function electron_plasma_frequency(density::NumberDensity)
    upreferred(sqrt(density * elementry_charge^2 / electron_mass / free_space_permittivity))
end

function ion_plasma_frequency(density::NumberDensity, charge_state::Integer, ion_mass::Unitful.Mass)
    upreferred(sqrt(density * (charge_state * elementry_charge)^2 / ion_mass / free_space_permittivity))
end

# electron and ion trapping rates excluded

#electron and ion collision rates excluded

function electron_debroglie_length(temp::Unitful.Temperature)
    upreferred(sqrt(2 * pi * reduced_planck_constant^2 / electron_mass / boltzmann_constant / temp))
end

function classical_minimum_approach_distance(temp::Unitful.Temperature)
    upreferred(elementry_charge^2 / boltzmann_constant / temp / (4 * pi * free_space_permittivity))
end

# TODO: electron_gyryradius

# TODO: ion_gyroradius

function electron_inertial_length(density::NumberDensity)
    upreferred(light_speed / electron_plasma_frequency(density))
end

function ion_inertial_length(density::NumberDensity, charge_state::Integer, ion_mass::Unitful.Mass)
    upreferred(light_speed / ion_plasma_frequency(density, charge_state, ion_mass))
end

function debye_length(density::NumberDensity, temp::Unitful.Temperature)
    upreferred(sqrt(free_space_permittivity * boltzmann_constant * temp / density / elementry_charge^2))
end

# TODO: does making energy and temperature interchangable make sense? Should I
# do this for everything?
function debye_length(density::NumberDensity, thermal_energy::Unitful.Energy)
    upreferred(sqrt(free_space_permittivity * thermal_energy / density / elementry_charge^2))
end

# TODO: plasma_parameter(temperature, density)

end
