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

    @test T96(param)(r_gsm, ps) ≈ [61.17831041891597, -1.4611958991749145, -40.44973158310239]
    @test T01(param)(r_gsm, ps) ≈ [46.972663449207076, 1.5442350206329172, -31.3541847716317] rtol = 1.0e-6
    @test TS04(param)(r_gsm, ps) ≈ [25.835474201385036, 1.5987724615979861, -18.1054945348421]

    # TA16_RBF.f (11/23/2021) compiled with gfortran; agrees to ~1e-9 as it truncates π and the degree-to-radian factor
    @test TA16(pdyn = 2.0, symhc = -30.0, xind = 0.8, byimf = 3.0)(r_gsm, ps) ≈ [36.193199414830254, 2.4386939542793806, -19.349898119812217] rtol = 1.0e-8
    @test TA16(pdyn = 6.0, symhc = -120.0, xind = 1.8, byimf = -7.0)((8.0, -3.0, 1.0), 0.4) ≈ GSM(54.775603850148237, -3.0340767211959898, 26.967360173238905) rtol = 1.0e-8

    @testset "near magnetopause and outside magnetosphere" begin
        # intermediate sigma case (S0 - DSIG < sigma < S0 + DSIG)
        @test t96([-6.5, 13, 13.0], ps, param...) ≈ [10.64621818721388, -0.8896267128450983, 1.9983159012100993]
        @test T01(param)([-6.5, 10, 13.0], ps) ≈ [33.90299388160385, -3.7496216888652825, -4.16584515175392] rtol = 1.0e-6
        @test T01(param)([-6.5, 10, 23.0], ps) ≈ [-0.281726292169127, 3.0187227770015053, -0.45328468009496703] rtol = 1e-14
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
    @test collect(TsyganenkoModels.shlcar(TsyganenkoModels.BIRK1_SHLD, Val(4), 1, 2, 3, 4)) ≈ [2.525027295079125, 0.05065465250867974, 0.12243989839075299]

    # Geopack.geopack.t96.birk1tot_02(1,2,3,-0.46)
    @test collect(TsyganenkoModels.birk1tot_02(1, 2, 3, -0.46)) ≈ [-4.138521384351386, 1.7274528180037634, -1.957000062050696]

    ps = -0.46049650108726486
    @test collect(TsyganenkoModels.dipole_shield(1, 2, 3, ps)) ≈ [-6.624367426893622, 0.5530405783442072, 19.261862546642284]

    @test collect(TsyganenkoModels.rc_symm(0, 0, 1)) ≈ [0.0, 0.0, -15.875017940770613] rtol = 1.0e-6
    @test collect(TsyganenkoModels.prc_quad(0, 0, 1)) ≈ [-38.33420080986639, 0.0, 0.0] rtol = 1.0e-6
end

@testset "Field line tracing" begin
    using Dates
    using GeoCotrans
    using OrdinaryDiffEqTsit5

    # Test case from https://github.com/tsssss/geopack/blob/master/geopack/test_geopack1.md
    # geopack dir=-1 = parallel to B (geopack sign convention is opposite to GeoCotrans dir=1)
    time = DateTime(2001, 1, 1, 2, 3, 4)
    r_gsm = GSM(-5.1, 0.3, 2.8)
    param = (; pdyn = 2.0, dst = -87.0, byimf = 2.0, bzimf = -5.0)

    # Geopack reference footpoints (parallel to B, r0=1.1 RE, w1-w6=0)
    # ~1% systematic offset vs geopack expected: geopack uses fixed-step RK4 (ds=0.5 RE, ~8 steps)
    # compared to our more accurate fully converged Tsit5
    footpoints = (;
        T89 = GSM(-0.7169607888154808, 0.03270815651368241, 0.8251160327006916),
        T96 = GSM(-0.7231220043155986, 0.03524984007351079, 0.818160665907448),
        T01 = GSM(-0.7206908533178673, 0.039122906017394675, 0.8206093486807777),
        TS04 = GSM(-0.7147665538572625, 0.035790564238260165, 0.8271336967350074),
    )

    models = (;
        T89 = TsyIGRF(T89(iopt = 2)),
        T96 = TsyIGRF(T96(param)),
        T01 = TsyIGRF(T01(param)),
        TS04 = TsyIGRF(TS04(param)),
    )

    for name in keys(models)
        sol = trace(r_gsm, time, Tsit5(); model = models[name], dir = 1, r0 = 1.1, rlim = 10.0)
        footpoint = GSM(sol.u[end]...)
        @test footpoint ≈ footpoints[name] rtol = 1.0e-2
    end
end

@testset "Driven models" begin
    using Dates
    using TsyganenkoModels: fill_gaps!, t01_g, ts05_w
    t = DateTime(2001, 1, 1, 2, 3, 4)
    r = GSM(-5.1, 0.3, 2.8)
    param = (; pdyn = 2.0, dst = -87.0, byimf = 2.0, bzimf = -5.0)
    @test TsyIGRF(Driven(T96, _ -> param))(r, t) == TsyIGRF(T96(param))(r, t)

    x = [NaN, 1.0, NaN, NaN, 4.0, NaN, NaN, NaN, 8.0, NaN]
    @test isequal(fill_gaps!(x; maxgap = 2), [NaN, 1, 2, 3, 4, NaN, NaN, NaN, 8, NaN])

    # Tsyganenko (2002b) §3.1: V = 400 km/s, By = 0, Bz = -5 nT gives G1 ≈ 6, G2 = 10
    n = 4000
    v, by, bz = fill(400.0, n), zeros(n), fill(-5.0, n)
    @test collect(t01_g(v, by, bz, n)) ≈ [400 * 0.125^2 / 1.125, 10]

    # Driving functions S_k = 1 (density 5/1.16, 400 km/s, Bs = 5 nT) held steady: W_k → 1
    den = fill(5 / 1.16, n)
    @test all(w -> isapprox(w, 1; rtol = 0.05), ts05_w(den, v, bz)[n])
    # Memory carries across a gap: a missing sample counts as no driving
    gappy, northward = copy(bz), copy(bz)
    gappy[n - 100] = NaN
    northward[n - 100] = 1
    w = ts05_w(den, v, gappy)
    @test w[n] == ts05_w(den, v, northward)[n]
    @test all(isnan, w[n - 100])

    # TA16 Newell index: By = -3, Bz = -4 nT has sin²(θ/2) = 0.9
    p = TsyganenkoModels.ta16_inputs(ones(40), zeros(40), fill(400.0, 40), fill(-3.0, 40), fill(-4.0, 40), 20)
    @test p.xind ≈ 1.0e-4 * 400^(4 / 3) * 5^(2 / 3) * 0.9^(4 / 3)
    # Trailing window: 5-min By of 2015-02-28 03:10–03:40 and ⟨By⟩ = -0.71 at 03:40 from Tsyganenko's
    # 2015_OMNI_5m_with_GEO_RBF_TA17_drivers.dat; 6-sample windows give -1.10 or -0.37.
    by = [1.58, 1.17, 1.17, -0.68, -2.53, -2.91, -2.79, 0, 0, 0]
    @test TsyganenkoModels.ta16_inputs(ones(10), zeros(10), fill(400.0, 10), by, zeros(10), 7).byimf ≈ -0.71 atol = 0.005
end

@testset "OMNI input windows" begin
    using CDAWeb, Dates
    using TsyganenkoModels: inputs
    ext = Base.get_extension(TsyganenkoModels, :TsyganenkoModelsCDAWebExt)
    day = Date(2001, 1, 1)  # synthetic data, not fetched
    try
        k = Float64.((1 - ext.PAD):(ext.PER_DAY + ext.PAD))
        n = length(k)
        pdyn = 100 .+ k
        pdyn[ext.PAD + 2] = NaN
        ext.HOURLY[(day, "KP")] = fill(30.0, 24)
        ext.HOURLY[(day, "KP")][2] = NaN
        ext.HOURLY[(day, "DST")] = collect(1.0:24)
        t = DateTime(day)
        @test inputs(Driven(T89), t).iopt === 4
        @test_throws ArgumentError inputs(Driven(T89), t + Hour(1))
        @test !haskey(ext.OMNI_DAYS, day)
        ext.OMNI_DAYS[day] = ext.OMNIDay(pdyn, 400 .+ k, zeros(n), zeros(n), fill(-5.0, n), k)
        ext.TS04_W[day] = [(x, x, x, x, x, x) for x in 1.0:ext.PER_DAY]
        ext.TA16_OBS[day] = ext.TA16Obs(zeros(n), zeros(n), fill((-400.0, 0.0, 0.0), n))
        for minutes in (59, 60, 1439)
            h = minutes ÷ 60 + 1
            p = inputs(Driven(T96), t + Minute(minutes))
            @test p.pdyn ≈ (h == 1 ? 100 + 76 / 11 : 100 + 12h - 5.5)
            @test p.dst == h
        end
        for minutes in (0, 1439)
            sample = minutes ÷ 5 + 1
            p = inputs(Driven(T01), t + Minute(minutes))
            @test (p.pdyn, p.dst) == (100 + sample, sample)
            @test p.g2 ≈ 0.025 * (400 + sample - 5.5)
            @test inputs(Driven(TS04), t + Minute(minutes)).w6 == sample
        end
        @test_throws ArgumentError inputs(Driven(T01), t + Minute(5))
        for minutes in (0, 1439)
            sample = minutes ÷ 5 + 1
            p = inputs(Driven(TA16), t + Minute(minutes))
            valid = filter(!=(2), collect((sample - 3):(sample + 3)))
            @test p.symhc ≈ sum(0.8j - 13sqrt(100 + j) for j in valid) / length(valid)
        end
    finally
        foreach(c -> delete!(c, day), (ext.OMNI_DAYS, ext.TS04_W, ext.TA16_OBS))
        foreach(name -> delete!(ext.HOURLY, (day, name)), ("KP", "DST"))
    end
end

@testset "OMNI-driven models" begin
    using CDAWeb, Dates
    using TsyganenkoModels: inputs
    t = DateTime(2017, 9, 8, 12, 2)  # Kp = 8- storm
    @test inputs(Driven(T89), t).iopt == 7
    r = (-5.1, 0.3, 2.8)
    for M in (T96, T01, TS04)
        m = Driven(M)
        @test m(r, t) == M(inputs(m, t))(r, dipole_tilt(t))
    end
    using GeoCotrans: GSW, at, rotation
    p = inputs(Driven(TA16), t)
    f = GSW(p.vgse)
    R = rotation(GSM() => f, t)
    B = R' * TA16(; p...)(R * collect(r), dipole_tilt(t, f))
    @test Driven(TA16)(r, t) ≈ B
    @test at(Driven(TA16), t)(collect(r)) ≈ B
    # Pdyn of Tsyganenko's 2015_OMNI_5m_with_GEO_RBF_TA17_drivers.dat at 2015-04-15 07:05
    @test inputs(Driven(TA16), DateTime(2015, 4, 15, 7, 7)).pdyn ≈ 8.27 atol = 0.005
    quiet(f) = Base.CoreLogging.with_logger(f, Base.CoreLogging.NullLogger())
    @test_throws ArgumentError quiet(() -> inputs(Driven(T96), DateTime(1960)))
end
