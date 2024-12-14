"""
    plasma_beta(T, n, B)

Compute the plamsa beta (β), the ratio of thermal pressure to magnetic pressure.

# Arguments
- `T`: The temperature of the plasma.
- `n`: The particle density of the plasma.
- `B`: The magnetic field in the plasma.
"""
function plasma_beta(T::EnergyOrTemp, n::NumberDensity, B::Unitful.BField)
    p_th = thermal_pressure(T, n)
    p_B = magnetic_pressure(B)
    return p_th / p_B |> upreferred
end