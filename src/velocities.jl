"""
    diamagnetic_drift(dp, B::BFields, n, q)

Calculate the diamagnetic drift velocity given by:

```math
v = - (∇p × B) / (q n |B|^2)
```

where ``∇p`` is the pressure gradient.
"""
function diamagnetic_drift(∇p, B::BFields, n, q = Unitful.q)
    -(∇p × B) / (q * n * B ⋅ B) .|> upreferred
end

"""
    ExB_drift(E, B)

Calculate the ``E × B`` drift velocity given by:

```math
v = (E × B) / |B|^2
```
"""
function ExB_drift end
@permute_args ExB_drift(E::EFields, B::BFields) = E × B / (B ⋅ B) .|> upreferred

"""
    force_drift(F, B::BFields, q)

Calculate the general force drift for a particle in a magnetic field given by:
```math
v = (F × B) / (q * |B|^2)
```
"""
force_drift(F, B::BFields, q = Unitful.q) = F × B / (q * B ⋅ B) .|> upreferred