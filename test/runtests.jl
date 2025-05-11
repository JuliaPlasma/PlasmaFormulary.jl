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
end

@testitem "dimensionless.jl" setup = [Py] begin
    using Unitful

    T = 2.0u"eV"
    n = 1e19u"m^-3"
    B = 0.1u"T"
    @test plasma_beta(T, n, B) ==
          plasma_beta(n, B, PlasmaFormulary.temperature(T)) ==
          plasma_beta(n, T, B)

    @test pyustrip(
        units.m,
        formulary.lengths.Debye_length(10 * units.eV, 1e19 / units.m^3),
    ) ≈ ustrip(u"m", PlasmaFormulary.debye_length(1e19 * u"m^-3", 10 * u"eV")) rtol = 1e-3
end
