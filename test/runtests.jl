using TestItems, TestItemRunner
using CondaPkg
CondaPkg.add(["astropy", "plasmapy"])

@testsnippet Py begin
    using PythonCall
    include("test_utils.jl")
    units = pyimport("astropy.units")
    formulary = pyimport("plasmapy.formulary")
end

@run_package_tests

@testitem "Misc" begin
    using Unitful
    using LinearAlgebra
    @test gyroradius(0.2u"T", Unitful.mp, Unitful.q, 1e6u"K") ≈ 0.0067u"m" rtol = 0.01
    @test_broken electron_gyroradius(1e6u"K", 0.2u"T", Unitful.q, Unitful.mp) ≈ 0.0067u"m" rtol =
        0.01

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
end

@testitem "dimensionless.jl" setup = [Py] begin
    using Unitful

    T = 2.0u"eV"
    n = 1e19u"m^-3"
    B = 0.1u"T"
    @test plasma_beta(T, n, B) ==
          plasma_beta(n, B, PlasmaFormularyZ.temperature(T)) ==
          plasma_beta(n, T, B)

    @test pyustrip(
        units.m,
        formulary.lengths.Debye_length(10 * units.eV, 1e19 / units.m^3),
    ) ≈ ustrip(u"m", PlasmaFormularyZ.debye_length(1e19 * u"m^-3", 10 * u"eV")) rtol = 1e-3
end
