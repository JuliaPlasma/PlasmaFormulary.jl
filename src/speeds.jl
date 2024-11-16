# TODO: Sound speed

const BFields = Union{BField,AbstractVector{<:BField}}

"""
    Alfven_velocity(B::BFields, ρ::Density)
    Alfven_velocity(B::BFields, n::NumberDensity; mass_numb = 1, z = 1)

Calculate the Alfvén velocity (v_A = B/√(μ₀ρ)) in a magnetized plasma.

# Arguments
- `B`: Magnetic field strength (scalar or vector)
- `ρ`: Mass density
- `n`: Number density
- `mass_numb=1`: Mass number of the ion species
- `z=1`: Ion charge number

# Examples
```julia
# Using mass density
B = 1.0u"T"
ρ = 1.0u"kg/m^3"
va = Alfven_velocity(B, ρ)

# Using number density
n = 1.0e19u"m^-3"
va = Alfven_velocity(B, n, mass_numb=1, z=1)

# Vector input
B = [1.0, 2.0, 3.0]u"T"
va = Alfven_velocity(B, n)
```

See also: [`Alfven_speed`](@ref)
"""
function Alfven_velocity end

@permutable_args function Alfven_velocity(B::BFields, ρ::Density)
    return B ./ sqrt(μ0 * ρ) .|> upreferred
end

@permutable_args function Alfven_velocity(
    B::BFields,
    n::NumberDensity;
    mass_numb = 1,
    z = 1,
)
    ρ = n * (mass_numb * mp + z * me)
    return Alfven_velocity(B, ρ)
end

"""
    Alfven_speed(B::BField, ρ::Density)
    Alfven_speed(B::BField, n::NumberDensity; mass_numb = 1, z = 1)

The typical propagation speed of magnetic disturbances in a quasineutral plasma.

References: [Wikipedia](https://en.wikipedia.org/wiki/Alfven_wave), [PlasmaPy API](https://docs.plasmapy.org/en/latest/api/plasmapy.formulary.speeds.Alfven_speed.html)
"""
function Alfven_speed end

@permutable_args Alfven_speed(B::BFields, n::NumberDensity; mass_numb = 1, z = 1) =
    norm(Alfven_velocity(B, n; mass_numb, z))
@permutable_args Alfven_speed(B::BFields, ρ::Density) = norm(Alfven_velocity(B, ρ))


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