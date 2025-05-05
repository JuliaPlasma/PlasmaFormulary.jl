using Chairmarks
using PlasmaFormulary
using Random
using Unitful

B = 0.1u"T"
Bs = rand(100) * B

@b gyrofrequency(B, :p)
@b gyrofrequency.(Bs, :p)
@b gyrofrequency.(:p, Bs)