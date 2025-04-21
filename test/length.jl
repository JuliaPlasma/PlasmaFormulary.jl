@testitem "Length" begin
    using Unitful
    @test inertial_length(5u"m^-3", :e) ≈ 2376534.75u"m"
    @test inertial_length(5u"m^-3", :"He+") ≈ 2.02985e8u"m" rtol = 0.0001
end