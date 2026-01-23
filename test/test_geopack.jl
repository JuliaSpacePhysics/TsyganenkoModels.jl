using Test

@testset "Geopack.jl" begin
    using Geopack
    using Geopack: geopack
    using Geopack: pyconvert
    ps = Geopack.recalc(100)
    @test ps == -0.46049650108726486
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
    geopack.t04.rh0 = geopack.t01.rh0
    geopack.t04.g = geopack.t01.g
    geopack.t04.d = geopack.t01.d
    geopack.t04.dxshift1 = geopack.t01.dxshift1
    geopack.t04.dxshift2 = geopack.t01.dxshift2
    geopack.t04.deltady = geopack.t01.deltady

    # compare birk_tot
    @test pyconvert(Vector{Float64}, geopack.t04.shlcar3x3(1, 2, 3, ps)) == [-6.624367426893622, 0.5530405783442072, 19.261862546642284]
    @test pyconvert(Bool, geopack.t04.birk_tot(0, ps, 1, 2, 3) == geopack.t01.birk_tot(0, ps, 1, 2, 3))
    @test pyconvert(Bool, geopack.t04.shlcar3x3(1, 2, 3, ps) == geopack.t01.shlcar3x3(1, 2, 3, ps))
    @test pyconvert(Bool, geopack.t04.full_rc(0, ps, 1, 2, 3) == geopack.t01.full_rc(0, ps, 1, 2, 3))
    @test pyconvert(Bool, geopack.t04.deformed(0, ps, 1, 2, 3) == geopack.t01.deformed(0, ps, 1, 2, 3))

    xgsm, ygsm, zgsm = (-5.1, 0.3, 2.8)
    ps = -0.533585131
    @test Geopack.t04([pdyn, dst, byimf, bzimf, 0, 0, 0, 0, 0, 0, 0], ps, xgsm, ygsm, zgsm) == [25.835474201385036, 1.5987724615979861, -18.1054945348421]
end
