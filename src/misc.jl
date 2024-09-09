"""
    thermal_pressure(T, n)

Calculate the thermal pressure for a Maxwellian distribution.

# Arguments
- `T`: The particle temperature or energy.
- `n`: The particle number density.
"""
function thermal_pressure(T::EnergyOrTemp, n::NumberDensity)
    return n * k * temperature(T) |> upreferred
end

const p_th = thermal_pressure

"""
    magnetic_pressure(B)

Calculate the magnetic pressure.
"""
function magnetic_pressure(B::BField)
    return (B^2) / (2 * μ0)
end