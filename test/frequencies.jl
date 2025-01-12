@testset "Frequency" begin

    @test gyrofrequency(0.01u"T", :p) ≈ 957883.3292211705u"s^-1"
    @test plasma_frequency(1e19u"m^-3", :p) ≈ 4163294534.0u"s^-1"
    @test plasma_frequency(1e19u"m^-3") ≈ 1.783986365e11u"s^-1"
end