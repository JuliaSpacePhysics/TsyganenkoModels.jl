# Reference: https://github.com/tsssss/geopack/blob/master/geopack/test_geopack1.md

using TsyganenkoModels
using Test
using Aqua
using Chairmarks

@testset "Code quality (Aqua.jl)" begin
    Aqua.test_all(TsyganenkoModels)
end

@testset "External magnetic fields" begin
    using Dates
    time = DateTime(2001, 1, 1, 2, 3, 4)
    r_gsm = (-5.1, 0.3, 2.8)
    ps = -0.533585131
    iopt = 2
    B_true = [20.77213175686351, -0.6465547428023687, -15.071641970338984]
    @test T89(iopt)(r_gsm, ps) == B_true
    @test T89(; iopt)(r_gsm, time) ≈ B_true rtol = 1.0e-4

    param = (; pdyn = 2.0, dst = -87.0, byimf = 2.0, bzimf = -5.0)

    result = T96(param)(r_gsm, ps)
    @test result ≈ [61.17831041891597, -1.4611958991749145, -40.44973158310239]
    # T01 model test
    result_t01 = T01(param)(r_gsm, ps)
    @test result_t01 ≈ [46.972663449207076, 1.5442350206329172, -31.3541847716317] rtol = 1.0e-6
    @test TS04(param)(r_gsm, ps) ≈ [25.835474201385036, 1.5987724615979861, -18.1054945348421]

    @testset "near magnetopause and outside magnetosphere" begin
        # intermediate sigma case (S0 - DSIG < sigma < S0 + DSIG)
        @test t96([-6.5, 13, 13.0], ps, param...) ≈ [10.64621818721388, -0.8896267128450983, 1.9983159012100993]
        @test T01(param)([-6.5, 10, 13.0], ps) ≈ [33.90299388160385, -3.7496216888652825, -4.16584515175392] rtol = 1.0e-6
        @test T01(param)([-6.5, 10, 23.0], ps) == [-0.281726292169127, 3.0187227770015053, -0.45328468009496703]
    end
end


@testset "Composite magnetic fields" begin
    using Dates
    time = DateTime(2001, 1, 1, 2, 3, 4)
    r_gsm = GSM(-5.1, 0.3, 2.8)
    model = TsyIGRF()
    B_true = [289.34552144342604, -20.162245780356457, -68.44100172746695]
    @test IGRF()(r_gsm, time) + T89(iopt = 3)(r_gsm, time) ≈ B_true
    @test model(r_gsm, time) ≈ B_true

    @info @b TsyIGRF()($r_gsm, $time), IGRF()($r_gsm, $time) + T89(iopt = 3)($r_gsm, $time)
end

@testset "External magnetic fields - Components" begin
    @test TsyganenkoModels.r2inner(1, 2, 3) == (-10.166350217481673, 16.22849226559209, -4.448363653306496)
    @test collect(TsyganenkoModels.r2sheet(1, 2, 3)) ≈ [9.453799383727832, 2.8367547754153284, 0.8914493552029839]

    # Geopack.geopack.t96.birk1shld(1,2,3,4)
    @test collect(TsyganenkoModels.birk1shld(1, 2, 3, 4)) ≈ [2.525027295079125, 0.05065465250867974, 0.12243989839075299]

    # Geopack.geopack.t96.birk1tot_02(1,2,3,-0.46)
    @test collect(TsyganenkoModels.birk1tot_02(1, 2, 3, -0.46)) ≈ [-4.138521384351386, 1.7274528180037634, -1.957000062050696]

    ps = -0.46049650108726486
    @test collect(TsyganenkoModels.shlcar3x3(1, 2, 3, ps)) ≈ [-6.624367426893622, 0.5530405783442072, 19.261862546642284]

    @test TsyganenkoModels.rc_symm(0, 0, 1) == (-0.0, -0.0, -15.875017940770613)
    @test TsyganenkoModels.prc_quad(0, 0, 1) == (-38.33420080986639, 0.0, 0.0)
end
