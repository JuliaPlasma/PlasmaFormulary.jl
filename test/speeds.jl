@testitem "Alfven speed & Alfven velocity" begin
    using LinearAlgebra: norm
    using Unitful

    B = [-0.014u"T", 0.0u"T", 0.0u"T"]
    Bmag = norm(B)
    n = 5e19 * u"m^-3"
    mass_number = 1
    rho = n * mass_number * Unitful.mp
    @test Alfven_speed(B, rho) ==
          Alfven_speed(rho, B) ==
          norm(Alfven_velocity(B, rho)) ==
          Alfven_velocity(Bmag, rho) ==
          Alfven_velocity(Bmag, n, mass_number) ≈
          43185.625u"m/s"

    @test Alfven_velocity(B, rho) ==
          Alfven_velocity(rho, B) ≈
          [-43185.625u"m/s", 0.0u"m/s", 0.0u"m/s"]
end

using Test
using Unitful
using PlasmaFormularyZ

@testitem "Ion sound speed" begin
    using Unitful
    using Unitful: mp
    # Typical values for hydrogen plasma
    T_e = 5e6u"K"
    T_i = 0u"K"
    m_i = mp  # proton mass
    Z = 1

    # Non-dispersive limit
    cs = ion_sound_speed(T_e, T_i, m_i, Z)
    @test isapprox(cs, 2.03155e5u"m/s"; rtol = 1e-5)

    # Dispersive limit, low k (should be close to non-dispersive)
    n_e = 5e19u"1/m^3"
    k1 = 3e1u"1/m"
    cs1 = ion_sound_speed(T_e, T_i, m_i, Z; n_e, k = k1)
    @test isapprox(cs1, 2.03155e5u"m/s"; rtol = 1e-5)

    # Dispersive limit, high k (should be much lower)
    k2 = 3e7u"1/m"
    cs2 = ion_sound_speed(T_e, T_i, m_i, Z; n_e, k = k2)
    @test cs2 < cs1

    # With nonzero T_i
    cs3 = ion_sound_speed(500u"eV", 200u"eV", 3.343583719e-27u"kg", 1; n_e = n_e, k = k1)
    @test isapprox(cs3, 2.29585e5u"m/s"; rtol = 1e-5)
end