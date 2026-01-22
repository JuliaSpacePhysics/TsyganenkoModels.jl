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
    @test t89(r_gsm, time, iopt) ≈ B_true rtol = 1e-4

    pdyn = 2.0   # Solar wind dynamic pressure [nPa]
    dst = -87.0  # Dst index [nT]
    byimf = 2.0  # IMF By [nT]
    bzimf = -5.0 # IMF Bz [nT]
    result = t96(r_gsm, ps, pdyn, dst, byimf, bzimf)
    @test collect(result) ≈ [61.17831041891597, -1.4611958991749145, -40.44973158310239]
end
