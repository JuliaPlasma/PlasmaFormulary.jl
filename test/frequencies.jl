@testitem "Frequency" begin
    using DimensionfulAngles
    using Unitful
    @test gyrofrequency(0.01u"T", :p) ≈ 957883.3292211705u"radᵃ*s^-1"
    @test uconvert(u"Hz", gyrofrequency(0.1u"T", :e), Periodic()) ≈ 2.799248987233304e9u"Hz"
    @test plasma_frequency(1e19u"m^-3", :p) ≈ 4163294534.0u"radᵃ*s^-1"
    @test plasma_frequency(1e19u"m^-3") ≈ 1.783986365e11u"radᵃ*s^-1"
end