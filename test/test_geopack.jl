using Test

@testset "Geopack.jl" begin
    using Geopack
    using Geopack: geopack
    using Geopack: pyconvert
    ps = Geopack.recalc(100)
    pdyn = 2.0; dst = -87.0
    byimf = 2.0; bzimf = -5.0
    g1 = 0.0; g2 = 0.0
    parmod = [pdyn, dst, byimf, bzimf, g1, g2]

    xgsm, ygsm, zgsm = 𝐫 = [1, 2, 3]
    @test Geopack.t01(parmod, ps, xgsm, ygsm, zgsm) == [20.302274127835044, -5.466320796811376, -66.92483279581937]
    # set model parameters
    geopack.t04.xkappa1 = geopack.t01.xkappa1
    geopack.t04.xkappa2 = geopack.t01.xkappa2
    geopack.t04.sc_sy = geopack.t01.sc_sy
    geopack.t04.sc_pr = geopack.t01.sc_pr
    geopack.t04.phi = geopack.t01.phi

    # compare birk_tot
    @test pyconvert(Bool, geopack.t04.birk_tot(0, ps, 1, 2, 3) == geopack.t01.birk_tot(0, ps, 1, 2, 3))
    @test pyconvert(Bool, geopack.t04.shlcar3x3(1, 2, 3, ps) == geopack.t01.shlcar3x3(1, 2, 3, ps))
    @test pyconvert(Bool, geopack.t04.full_rc(0, ps, 1, 2, 3) == geopack.t01.full_rc(0, ps, 1, 2, 3))
end