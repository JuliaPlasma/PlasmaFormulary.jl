@testitem "Current density" begin
    using Unitful
    using PlasmaFormulary: CurrentDensity
    @test Alfven_current_density(1u"T", 1u"m^-3") ≈ 0.0034946731981915545 * u"A/m^2"
    @test Alfven_current_density(1u"T", 1u"m^-3") isa CurrentDensity
    @test Alfven_current_density(1u"T", 1u"m") isa CurrentDensity
end