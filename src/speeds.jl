# TODO: Sound speed

"""
The typical propagation speed of magnetic disturbances in a quasineutral plasma.

References: [PlasmaPy API Documentation](https://docs.plasmapy.org/en/stable/api/plasmapy.formulary.speeds.Alfven_speed.html)
"""
function alfven_velocity(B::BField, ρ::Density)
    return B / sqrt(μ0 * ρ) |> upreferred
end

function alfven_velocity(B::BField, n::NumberDensity, mass_number=1)
    ρ = n * mass_number * mp
    return alfven_velocity(B, ρ)
end


"""
Calculate the Alfven velocity for magnetic field vector.
"""
alfven_velocity(𝐁::AbstractVector{<:BField}, ρ::Density) =
    map(𝐁) do B
        alfven_velocity(B, ρ)
    end

function alfven_velocity(𝐁::AbstractVector{<:BField}, n::NumberDensity, mass_number=1)
    alfven_velocity(𝐁, n * mass_number * mp)
end

# TODO: Add docstrings
abstract type ThermalVelocityMethod end
struct MostProbable <: ThermalVelocityMethod end
struct RMS <: ThermalVelocityMethod end
struct MeanMagnitude <: ThermalVelocityMethod end
struct NRL <: ThermalVelocityMethod end

"""
    thermal_velocity_coefficients(method::ThermalVelocityMethod, ndim::Int)

Get the thermal speed coefficient corresponding to the desired thermal speed definition.

# Arguments
- `method::ThermalVelocityMethod`: Method to be used for calculating the thermal speed. Valid values are `MostProbable()`, `RMS()`, `MeanMagnitude()`, and `NRL()`.
- `ndim::Val{Int}`: Dimensionality (1D, 2D, 3D) of space in which to calculate thermal speed. Valid values are `Val(1)`, `Val(2)`, or `Val{3}`.
"""
thermal_velocity_coefficients(::MostProbable, ::Val{1}) = 0.0
thermal_velocity_coefficients(::MostProbable, ::Val{2}) = 1.0
thermal_velocity_coefficients(::MostProbable, ::Val{3}) = sqrt(2)
thermal_velocity_coefficients(::RMS, ::Val{1}) = 1.0
thermal_velocity_coefficients(::RMS, ::Val{2}) = sqrt(2)
thermal_velocity_coefficients(::RMS, ::Val{3}) = sqrt(3)
thermal_velocity_coefficients(::MeanMagnitude, ::Val{1}) = sqrt(2 / π)
thermal_velocity_coefficients(::MeanMagnitude, ::Val{2}) = sqrt(π / 2)
thermal_velocity_coefficients(::MeanMagnitude, ::Val{3}) = sqrt(8 / π)
thermal_velocity_coefficients(::NRL, ::Val{1}) = 1.0
thermal_velocity_coefficients(::NRL, ::Val{2}) = 1.0
thermal_velocity_coefficients(::NRL, ::Val{3}) = 1.0

function thermal_velocity(
    T::EnergyOrTemp,
    mass::Unitful.Mass,
    method::ThermalVelocityMethod=MostProbable(),
    ndim=3,
)
    coeff = thermal_velocity_coefficients(method, Val(ndim))
    return coeff * sqrt(k * temperature(T) / mass)
end

function thermal_temperature(
    V::Unitful.Velocity,
    mass::Unitful.Mass,
    method::ThermalVelocityMethod=MostProbable(),
    ndim=3,
)
    coeff = thermal_velocity_coefficients(method, Val(ndim))
    return mass * V^2 / (k * coeff^2) |> upreferred
end

function electron_thermal_velocity(eot::EnergyOrTemp)
    thermal_velocity(eot, me, NRL())
end
