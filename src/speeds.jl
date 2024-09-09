# TODO: Sound speed

function Alfven_velocity(B::BField, ρ::Density)
    return B / sqrt(μ0 * ρ) |> upreferred
end

function Alfven_velocity(B::BField, n::NumberDensity; mass_numb = 1, Z = 1)
    ρ = n * (mass_numb * mp + Z * me)
    return Alfven_velocity(B, ρ)
end

Alfven_velocity(B::AbstractVector, args...; kwargs...) =
    [Alfven_velocity(Bi, args...; kwargs...) for Bi in B]

"""
    Alfven_speed(args...; kwargs...)

The typical propagation speed of magnetic disturbances in a quasineutral plasma.

References: [PlasmaPy API Documentation](https://docs.plasmapy.org/en/stable/api/plasmapy.formulary.speeds.Alfven_speed.html)
"""
Alfven_speed(args...; kwargs...) = norm(Alfven_velocity(args...; kwargs...))


function thermal_speed(
    T::EnergyOrTemp,
    mass::Unitful.Mass,
    method = "most_probable",
    ndim = 3,
)
    coeff = thermal_speed_coefficients(method, ndim)
    return coeff * sqrt(k * temperature(T) / mass)
end

thermal_speed(T::EnergyOrTemp, mass_numb, args...) =
    thermal_speed(T, mass_numb * u, args...)

function thermal_temperature(
    V::Unitful.Velocity,
    mass::Unitful.Mass,
    method = "most_probable",
    ndim = 3,
)
    coeff = thermal_speed_coefficients(method, ndim)
    return mass * V^2 / (k * coeff^2) |> upreferred
end

"""
    thermal_speed_coefficients(method::String, ndim::Int)

Get the thermal speed coefficient corresponding to the desired thermal speed definition.

# Arguments
- `method::String`: Method to be used for calculating the thermal speed. Valid values are `"most_probable"`, `"rms"`, `"mean_magnitude"`, and `"nrl"`.
- `ndim::Int`: Dimensionality (1D, 2D, 3D) of space in which to calculate thermal speed. Valid values are `1`, `2`, or `3`.
"""
function thermal_speed_coefficients(method::String, ndim::Int)
    # Define the coefficient dictionary
    _coefficients = Dict(
        (1, "most_probable") => 0.0,
        (2, "most_probable") => 1.0,
        (3, "most_probable") => sqrt(2),
        (1, "rms") => 1.0,
        (2, "rms") => sqrt(2),
        (3, "rms") => sqrt(3),
        (1, "mean_magnitude") => sqrt(2 / π),
        (2, "mean_magnitude") => sqrt(π / 2),
        (3, "mean_magnitude") => sqrt(8 / π),
        (1, "nrl") => 1.0,
        (2, "nrl") => 1.0,
        (3, "nrl") => 1.0,
    )

    # Attempt to retrieve the coefficient
    try
        return _coefficients[(ndim, method)]
    catch
        throw(
            ArgumentError(
                "Value for (ndim, method) pair not valid, got '($ndim, $method)'.",
            ),
        )
    end
end


function electron_thermal_velocity(eot::EnergyOrTemp)
    upreferred(sqrt(k * temperature(eot) / me))
end

function ion_thermal_velocity(eot::EnergyOrTemp, ion_mass::Unitful.Mass)
    upreferred(sqrt(k * temperature(eot) / ion_mass))
end