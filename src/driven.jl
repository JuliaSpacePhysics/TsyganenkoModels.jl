"""
    Driven(M, source = OMNI())

Model `M` evaluated at time `t` as `M(inputs(m, t))`.
"""
struct Driven{M, S} <: ExternalFieldModel
    source::S
end

Driven(M::Type{<:TsyganenkoModel}, source = OMNI()) = Driven{M, typeof(source)}(source)

native_csys(::Driven) = (GSM(), Cartesian3())
GeoCotrans.freeze(m::Driven{M}, t) where {M} = GeoCotrans.freeze(M(inputs(m, t)), t)

inputs(m::Driven{M}, t) where {M} = inputs(m.source, M, t)
inputs(f, ::Type, t) = f(t)

# The model's frame may depend on its inputs at `t` (TA16's GSW)
function evaluate_model(m::Driven, pos, t, from, ::Any, to; kw...)
    model = GeoCotrans.freeze(m, t)
    return evaluate_model(model, pos, t, from, native_csys(model), to; kw...)
end

"""
    OMNI()

Parameters from 5-min (and, for Kp and Dst, hourly) OMNI data:

- `pdyn = ρV²` with `ρ = 1.16 mₚN` (He/H = 0.04); gaps of at most 3 hours linearly interpolated, as Tsyganenko prepares TS05 inputs.
- `T89`: `iopt` from Kp. `T01`, `TS04`: SYM-H as `dst`, the index their fits used.
- `T96`: Dst and the means of the valid `pdyn`, `byimf`, `bzimf` samples over the UT hour containing `t`, as its fit tagged each field average with the concurrent hourly OMNI values ([Tsyganenko, 1996](https://ui.adsabs.harvard.edu/abs/1996ESASP.389..181T); [Tsyganenko, 2002b, ¶6](https://doi.org/10.1029/2001JA000220)).
- `T01`: other inputs at the 5-min sample containing `t`, the resolution of its fit, and `g1`, `g2` averaged over the valid samples in the hour ending there.
- `TS04`: `w1`–`w6` accumulated from 13 days before `t`, when the slowest memory (W3) has decayed below `e⁻¹⁰`; samples in unfilled gaps contribute nothing. Tsyganenko's yearly files instead restart the sums after each gap at the first 2-hour quiet period, and on January 1, discarding the decayed memory of earlier driving that the exponential kernel models.
- `TA16`: `pdyn` at the 5-min sample containing `t`; `symhc`, the mean of the valid `0.8 SYM-H - 13√pdyn` samples within ±15 min; `xind` and `byimf`, the means of the valid Newell N-index and By samples in the 30 min ending there ([Tsyganenko & Andreeva, 2016, eq. 9 and §4](https://doi.org/10.1002/2016JA023217)); `vgse`, the velocity at that sample with Earth's orbital 29.78 km/s added back to OMNI's Vy, defining the model's GSW frame, in which the IMF is taken. Tsyganenko's driver files use OMNI's Vy as is, which fits later THEMIS and Van Allen Probes magnetometer data worse.
"""
struct OMNI end

inputs(::OMNI, ::Type, t) = error("OMNI inputs require `using CDAWeb`")

# Linearly interpolate runs of at most `maxgap` NaNs bounded by valid values on both sides
function fill_gaps!(x; maxgap = 35)
    i = firstindex(x)
    while true
        j = findnext(isnan, x, i)
        isnothing(j) && return x
        k = findnext(!isnan, x, j)
        isnothing(k) && return x
        if j > firstindex(x) && k - j <= maxgap
            a, b = x[j - 1], x[k]
            for l in j:(k - 1)
                x[l] = a + (b - a) * (l - j + 1) / (k - j + 1)
            end
        end
        i = k
    end
end

function masked_mean(itr)
    s, n = 0.0, 0
    for x in itr
        isnan(x) && continue
        s += x
        n += 1
    end
    return s / n
end

# TA16 (symhc, xind, byimf) at sample `i` of 5-min series of `pdyn`, `symh`, solar wind speed `v` and IMF `by`, `bz` in the model frame.
# 30-min windows span 7 samples, both ends included:
# Tsyganenko's ⟨By⟩ in his archived TA17 driver files
# (https://web.archive.org/web/20230922174724/https://geo.phys.spbu.ru/~tsyganenko/models/geo/2015_OMNI_5m_with_GEO_RBF_TA17_drivers.dat)
# is the mean of By at t - 30 min, …, t.
function ta16_inputs(pdyn, symh, v, by, bz, i)
    trailing = (i - 6):i
    newell_n(j) = 1.0e-4 * v[j]^(4 / 3) * hypot(by[j], bz[j])^(2 / 3) * abs(sin(atan(by[j], bz[j]) / 2))^(8 / 3)
    return (;
        symhc = masked_mean(0.8 * symh[j] - 13 * sqrt(pdyn[j]) for j in (i - 3):(i + 3)),
        xind = masked_mean(newell_n(j) for j in trailing),
        byimf = masked_mean(view(by, trailing)),
    )
end

# T01 (G1, G2) at sample `i` of 5-min series of solar wind speed `v` [km/s] and IMF `by`, `bz` (GSM) [nT].
# T01 eq. (1), (4), Tsyganenko (2002b), https://doi.org/10.1029/2001JA000220: 12 preceding 5-min samples; we include the current one,
# as the TS05 W sums do.
function t01_g(v, by, bz, i)
    g1 = g2 = 0.0
    nvalid = 0
    for j in (i - 11):i
        (isnan(v[j]) || isnan(by[j]) || isnan(bz[j])) && continue
        nvalid += 1
        bperp = hypot(by[j], bz[j])
        h = (bperp / 40)^2 / (1 + bperp / 40)
        sin3 = bperp == 0 ? 0.0 : ((1 - bz[j] / bperp) / 2)^1.5  # sin³(θ/2), θ the IMF clock angle
        g1 += v[j] * h * sin3
        g2 += 0.005 * v[j] * max(-bz[j], 0.0)
    end
    return (g1 / nvalid, g2 / nvalid)
end

# Prepare_input_4.for and Parameters.par (A(45)–A(68)), Tsyganenko & Sitnov (2005), from
# https://geo.phys.spbu.ru/~tsyganenko/TS05_data_and_stuff/ (unreachable 2026-10; archived copies:
# https://web.archive.org/web/20090505024458/http://geo.phys.spbu.ru/~tsyganenko/TS05_data_and_stuff/Prepare_input_4.for
# https://web.archive.org/web/20220630063551/http://geo.phys.spbu.ru/~tsyganenko/TS05_data_and_stuff/Parameters.par)
const TS05_R = (0.383403, 0.648176, 0.0318752, 0.581168, 1.15070, 0.843004)  # relaxation rates [1/hour]
const TS05_EXP = (  # exponents of N, V, Bs
    (0.394732, 0.846509, 0.916555),
    (0.550920, 0.180725, 0.898772),
    (0.387365, 2.26596, 1.29123),
    (0.436819, 1.28211, 1.33199),
    (0.405553, 1.62290, 0.699074),
    (1.26131, 2.42297, 0.537116),
)
# TS05 (W1, …, W6) at each sample of 5-min series of proton density `n` [cm⁻³], solar wind speed `v` [km/s] and IMF `bz` (GSM) [nT],
# from zero memory at the first sample; NaN samples contribute nothing and get NaN.
function ts05_w(n, v, bz)
    δ = map(r -> r * 5 / 60, TS05_R)
    w = ntuple(_ -> 0.0, Val(6))
    out = Vector{NTuple{6, Float64}}(undef, length(n))
    for i in eachindex(n, v, bz, out)
        valid = !(isnan(n[i]) || isnan(v[i]) || isnan(bz[i]))
        w = map(w, δ, TS05_EXP) do wk, δk, (λ, β, γ)
            # 1.16: helium contribution, assuming He/H = 0.04
            s = valid && bz[i] < 0 ? (1.16 * n[i] / 5)^λ * (v[i] / 400)^β * (-bz[i] / 5)^γ : 0.0
            exp(-δk) * wk + δk * s
        end
        out[i] = valid ? w : ntuple(_ -> NaN, Val(6))
    end
    return out
end
