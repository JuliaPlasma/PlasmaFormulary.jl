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