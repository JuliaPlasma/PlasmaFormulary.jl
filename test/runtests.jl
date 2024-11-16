using PlasmaFormulary
using Test
using Aqua

@testset "PlasmaFormulary.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(PlasmaFormulary)
    end

    @testset "Calculations" begin
        include("test.jl")
    end

    @testset "Utils" begin
        include("utils.jl")
    end
end
