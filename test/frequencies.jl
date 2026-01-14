@testitem "Frequency" begin
    using DimensionfulAngles
    using Unitful
    @test plasma_frequency(1e19u"m^-3", :p) ≈ 4163294534.0u"radᵃ*s^-1"
    @test plasma_frequency(1e19u"m^-3") ≈ 1.783986365e11u"radᵃ*s^-1"
end