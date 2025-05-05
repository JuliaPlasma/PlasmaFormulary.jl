"""
    diamagnetic_drift(∇p, B::BFields, n, q)

Calculate the diamagnetic drift velocity given by:

```math
𝐯 = - (∇p × 𝐁) / (q n |𝐁|^2)
```

where ``∇p`` is the pressure gradient.
"""
function diamagnetic_drift(∇p, B::BFields, n, q = Unitful.q)
    -(∇p × B) / (q * n * B ⋅ B) .|> upreferred
end

"""
    ExB_drift(𝐄::EFields, 𝐁::BFields)

Calculate the ``E × 𝐁`` drift velocity given by:

```math
𝐯 = (𝐄 × 𝐁) / |𝐁|^2
```
"""
function ExB_drift end
@permute_args ExB_drift(𝐄::EFields, 𝐁::BFields) = 𝐄 × 𝐁 / (𝐁 ⋅ 𝐁) .|> upreferred

"""
    force_drift(𝐅, 𝐁::BFields, q)

Calculate the general force drift for a particle in a magnetic field given by:
```math
𝐯 = (𝐅 × 𝐁) / (q |𝐁|^2)
```
"""
force_drift(𝐅, 𝐁::BFields, q = Unitful.q) = 𝐅 × 𝐁 / (q * 𝐁 ⋅ 𝐁) .|> upreferred