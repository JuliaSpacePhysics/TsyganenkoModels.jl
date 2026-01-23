# Reference: https://github.com/tsssss/geopack/blob/master/geopack/test_geopack1.md

using TsyganenkoModels
using Test
using Aqua

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
    @test t89(r_gsm, ps, iopt) == B_true
    @test t89(r_gsm, time, iopt) ≈ B_true rtol = 1.0e-4

    pdyn = 2.0   # Solar wind dynamic pressure [nPa]
    dst = -87.0  # Dst index [nT]
    byimf = 2.0  # IMF By [nT]
    bzimf = -5.0 # IMF Bz [nT]
    result = t96(r_gsm, ps, pdyn, dst, byimf, bzimf)
    @test collect(result) ≈ [61.17831041891597, -1.4611958991749145, -40.44973158310239]

    # T01 model test (reference from Python geopack.t01)
    result_t01 = t01(r_gsm, ps, pdyn, dst, byimf, bzimf)
    @test collect(result_t01) ≈ [46.972663449207076, 1.5442350206329172, -31.3541847716317] rtol = 1.0e-6
    @test collect(ts04(r_gsm, ps, (pdyn, dst, byimf, bzimf))) ≈ [25.835474201385036, 1.5987724615979861, -18.1054945348421]
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
end
