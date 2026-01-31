using Test
using Geopack.PythonCall

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

    @test pyconvert(Tuple, geopack.t01.rc_symm(0, 0, 1)) == (-0.0, -0.0, -15.875017940770613)
    @test pyconvert(Tuple, geopack.t01.prc_quad(0, 0, 1)) == (-38.33420080986639, 0.0, 0.0)

    xgsm, ygsm, zgsm = r_gsm = [-5.1, 0.3, 2.8]
    ps = -0.533585131
    @test Geopack.t04([pdyn, dst, byimf, bzimf, 0, 0, 0, 0, 0, 0, 0], ps, xgsm, ygsm, zgsm) == [25.835474201385036, 1.5987724615979861, -18.1054945348421]

    # Test intermediate sigma case (S0 - DSIG < sigma < S0 + DSIG)
    r = (-6.5, 13, 13.0)
    ps = -0.533585131
    param = (; pdyn = 2.0, dst = -87.0, byimf = 2.0, bzimf = -5.0)
    @test Geopack.t96(collect(param), ps, r...) == [10.64621818721388, -0.8896267128450983, 1.9983159012100993]
    @test Geopack.t01([param..., 0.0, 0.0], ps, -6.5, 10, 13.0) == [33.90299388160385, -3.7496216888652825, -4.16584515175392]
    @test Geopack.t01([param..., 0.0, 0.0], ps, -6.5, 10, 23.0) == [-0.281726292169127, 3.0187227770015053, -0.45328468009496703]

    # Field line tracing: geopack dir=-1 = parallel to B (opposite sign to GeoCotrans dir=1)
    # recalc(time) must precede trace to set IGRF coefficients and dipole tilt in geopack globals
    # parmod for T96/T01/T04: [pdyn, dst, byimf, bzimf, w1..w6=0] (10 elements)
    # T89 parmod is just iopt (integer)
    @testset "dipole_tilt" begin
        using Dates, TsyganenkoModels
        time = DateTime("2015-10-16T00:00:00")
        ps = TsyganenkoModels.dipole_tilt(time)
        @test Geopack.recalc(time) ≈ ps rtol = 1.0e-4
    end

    @testset "track" begin
        using Dates
        trace_time = DateTime(2001, 1, 1, 2, 3, 4)
        Geopack.recalc(trace_time)
        parmod = [2., -87., 2., -5., 0., 0., 0., 0., 0., 0.]
        fp_t89  = pyconvert(Vector, geopack.geopack.trace(r_gsm..., -1, 10., 1.1, 2,      "t89",  "igrf"))
        fp_t96  = pyconvert(Vector, geopack.geopack.trace(r_gsm..., -1, 10., 1.1, parmod, "t96",  "igrf"))
        fp_t01  = pyconvert(Vector, geopack.geopack.trace(r_gsm..., -1, 10., 1.1, parmod, "t01",  "igrf"))
        fp_t04  = pyconvert(Vector, geopack.geopack.trace(r_gsm..., -1, 10., 1.1, parmod, "t04",  "igrf"))
        @test fp_t89[1:3] ≈ [-0.7169607888154808, 0.03270815651368241,  0.8251160327006916] rtol = 1e-6
        @test fp_t96[1:3] ≈ [-0.7231220043155986, 0.03524984007351079,  0.818160665907448]  rtol = 1e-6
        @test fp_t01[1:3] ≈ [-0.7206908533178673, 0.039122906017394675, 0.8206093486807777] rtol = 1e-6
        @test fp_t04[1:3] ≈ [-0.7147665538572625, 0.035790564238260165, 0.8271336967350074] rtol = 1e-6
    end
end
