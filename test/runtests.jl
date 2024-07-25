using PlasmaFormulary
using Test
using Aqua

@testset "PlasmaFormulary.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(PlasmaFormulary)
    end
    # Write your tests here.
end
