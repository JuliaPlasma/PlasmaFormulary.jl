using PlasmaFormulary
using Test
using Aqua
using Unitful

@testset "PlasmaFormulary.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(PlasmaFormulary)
    end

    @testset "Conversions" begin
        @test isapprox(uconvert(u"erg/K", PlasmaFormulary.boltzmann_constant), 1.38065e-16 * u"erg/K")
    end
end
