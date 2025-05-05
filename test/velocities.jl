@testitem "Drift velocities" begin
    using Unitful
    using Unitful: me

    ex = [1, 0, 0]
    ey = [0, 1, 0]
    ez = [0, 0, 1]

    # Diamagnetic drift
    dp = [0, -1, 0]u"Pa/m"
    B = [0, 0, 1]u"T"
    n = 1u"1/m^3"
    q = 1u"C"
    @test diamagnetic_drift(dp, B, n, q) ≈
          diamagnetic_drift(B, dp, n, q) ≈
          [1, -0, -0]u"m/s"

    # ExB drift
    E = ex * u"V/m"
    B = ey * u"T"
    @test ExB_drift(E, B) ≈ ExB_drift(B, E) ≈ [0, 0, 1]u"m/s"
    @test ExB_drift(ex * u"V/m", ex * u"T") ≈ [0, 0, 0]u"m/s"
    @test ExB_drift(ex * u"V/m", 100ey * u"T") ≈ [0, 0, 0.01]u"m/s"

    # Force drift
    g0 = 9.80665u"m/s^2"
    F = -ez * g0 * me
    B = ex * 0.01u"T"
    @test force_drift(F, B) ≈ force_drift(B, F) ≈ [0.0, -5.5756984e-9, 0.0]u"m/s"
    @test force_drift(F, ez * 0.01u"T") ≈ [0.0, -0.0, 0.0]u"m/s"
    @test force_drift(F, ex * u"T") ≈ [0.0, -5.5756984e-11, 0.0]u"m/s"
end