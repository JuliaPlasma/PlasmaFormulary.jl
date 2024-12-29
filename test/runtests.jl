using PlasmaFormulary
using Test
using Aqua
using Unitful
using LinearAlgebra

@testset "PlasmaFormulary.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(PlasmaFormulary)
    end

    @testset "Misc" begin
        @test gyroradius(0.2u"T", Unitful.mp, Unitful.q, 1e6u"K") ≈ 0.0067u"m" rtol=0.01
        @test_broken electron_gyroradius(1e6u"K", 0.2u"T", Unitful.q, Unitful.mp) ≈ 0.0067u"m" rtol=0.01

        B = [-0.014u"T", 0.0u"T", 0.0u"T"]
        Bmag = norm(B)
        n = 5e19 * u"m^-3"
        mass_number = 1
        rho = n * mass_number * Unitful.mp
        @test norm(alfven_velocity.(B, rho)) ==
              alfven_velocity(Bmag, rho) ==
              alfven_velocity(Bmag, n, mass_number) ≈
              43185.625u"m/s"

        @test alfven_velocity.(B, rho) ≈ [-43185.625u"m/s", 0.0u"m/s", 0.0u"m/s"]
        @test plasma_frequency(1e19u"m^-3", Unitful.q, Unitful.mp) ≈ 4163294534.0u"s^-1"
        @test plasma_frequency(1e19u"m^-3", 1, Unitful.mp / Unitful.u) ≈ 4163294534.0u"s^-1"
        @test plasma_frequency(1e19u"m^-3") ≈ 1.783986365e11u"s^-1"
    end

    @testset "dimensionless.jl" begin
        T = 2.0u"eV"
        n = 1e19u"m^-3"
        B = 0.1u"T"
        @test plasma_beta(T, n, B) == plasma_beta(n, B, PlasmaFormulary.temperature(T)) == plasma_beta(n, T, B)

        units = pyimport("astropy.units")
        formulary = pyimport("plasmapy.formulary")
    end
end
