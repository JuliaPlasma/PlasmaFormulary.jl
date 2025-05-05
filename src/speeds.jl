# TODO: Sound speed

"""
    Alfven_speed(B::BField, ρ)
    Alfven_speed(𝐁::Vector{BField}, ρ)
    Alfven_speed(B::BField, n::NumberDensity, mass_number = 1)
    Alfven_speed(𝐁::Vector{BField}, n::NumberDensity, mass_number = 1)

Alfvén speed ``V_A``, the typical propagation speed of magnetic disturbances in a quasineutral plasma.

Note that this is different from the Alfven velocity, see also [`Alfven_velocity`](@ref).
References: [PlasmaPy API Documentation](https://docs.plasmapy.org/en/stable/api/plasmapy.formulary.speeds.Alfven_speed.html)
"""
function Alfven_speed end
@permute_args Alfven_speed(𝐁::Union{BField,BFields}, ρ::Density) =
    norm(𝐁) / sqrt(μ0 * ρ) |> upreferred
@permute_args Alfven_speed(𝐁::Union{BField,BFields}, n::NumberDensity, mass_number) =
    Alfven_speed(𝐁, n * mass_number * mp)
@permute_args Alfven_speed(𝐁::Union{BField,BFields}, n::NumberDensity; mass_number = 1) =
    Alfven_speed(𝐁, n, mass_number)

"""
    Alfven_velocity(B::BField, ρ)
    Alfven_velocity(𝐁::Vector{BField}, ρ)

Calculate the Alfven velocity for magnetic field vector. See also [`Alfven_speed`](@ref).
"""
function Alfven_velocity end
@permute_args Alfven_velocity(B::Union{BField,BFields}, ρ::Density) =
    @. B / sqrt(μ0 * ρ) |> upreferred
@permute_args Alfven_velocity(B::Union{BField,BFields}, n::NumberDensity, mass_number) =
    Alfven_velocity(B, n * mass_number * mp)
@permute_args Alfven_velocity(B::Union{BField,BFields}, n::NumberDensity; mass_number = 1) =
    Alfven_velocity(B, n, mass_number)


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
    method::ThermalVelocityMethod = MostProbable(),
    ndim = 3,
)
    coeff = thermal_velocity_coefficients(method, Val(ndim))
    return coeff * sqrt(k * temperature(T) / mass)
end

function thermal_temperature(
    V::Unitful.Velocity,
    mass::Unitful.Mass,
    method::ThermalVelocityMethod = MostProbable(),
    ndim = 3,
)
    coeff = thermal_velocity_coefficients(method, Val(ndim))
    return mass * V^2 / (k * coeff^2) |> upreferred
end

function electron_thermal_velocity(eot::EnergyOrTemp)
    thermal_velocity(eot, me, NRL())
end
