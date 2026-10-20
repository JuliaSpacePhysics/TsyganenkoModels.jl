module TsyganenkoModelsCDAWebExt

using CDAWeb: CDAWeb, getdata, times
using Dates
using TsyganenkoModels: TsyganenkoModels, OMNI, T89, T96, T01, TS04, fill_gaps!, hourly_mean, t01_g, ts05_w, TS05_R

const STEP = Minute(5)
const PER_DAY = 288
# 5-min samples fetched before the day, for W memory to decay below e⁻¹⁰ at the slowest rate (W3, ~13 days),
# and after it, for filling gaps that straddle midnight
const NBEFORE = ceil(Int, 10 / (minimum(TS05_R) * 5 / 60))
const NAFTER = 36

struct OMNIDay
    pdyn::Vector{Float64}
    symh::Vector{Float64}
    by::Vector{Float64}
    bz::Vector{Float64}
    sw_hourly::Vector{NTuple{3, Float64}}  # pdyn, by, bz per UT hour
    g::Vector{NTuple{2, Float64}}
    w::Vector{NTuple{6, Float64}}
    kp::Vector{Float64}
    dst::Vector{Float64}
end

# Values of `name` on the regular grid `t0 + (0:n-1) * step`, NaN where OMNI has no record
function _series(ds, name, t0, n, step)
    x = getdata(ds[name], t0, t0 + n * step)
    return _grid!(fill(NaN, n), times(x), Array(x), t0, step)
end

# Function barrier: `getdata`'s return type is not inferred
function _grid!(out, ts, vs, t0, step)
    for (t, v) in zip(ts, vs)
        out[(DateTime(t) - t0) ÷ step + 1] = v
    end
    return out
end

function OMNIDay(day::Date)
    t0 = DateTime(day) - NBEFORE * STEP
    n = NBEFORE + PER_DAY + NAFTER
    hro = CDAWeb.Dataset("OMNI_HRO_5MIN")
    den, v, by, bz, symh = map(("proton_density", "flow_speed", "BY_GSM", "BZ_GSM", "SYM_H")) do name
        fill_gaps!(_series(hro, name, t0, n, STEP))
    end
    pdyn = @. 1.937e-6 * den * v^2  # 1.16 mₚ in nPa per (cm⁻³ (km/s)²), Tsyganenko's value
    idx = NBEFORE .+ (1:PER_DAY)
    hourly = CDAWeb.Dataset("OMNI2_H0_MRG1HR")
    return OMNIDay(
        pdyn[idx], symh[idx], by[idx], bz[idx],
        [map(x -> hourly_mean(x, NBEFORE + 12h), (pdyn, by, bz)) for h in 1:24],
        [t01_g(v, by, bz, i) for i in idx],
        ts05_w(den, v, bz)[idx],
        _series(hourly, "KP", DateTime(day), 24, Hour(1)),
        _series(hourly, "DST", DateTime(day), 24, Hour(1)),
    )
end

const CACHE = Dict{Date, OMNIDay}()
const LOCK = ReentrantLock()

function _sample(t)
    t = DateTime(t)
    day = Date(t)
    d = @lock LOCK get!(() -> OMNIDay(day), CACHE, day)
    return d, (t - DateTime(day)) ÷ STEP + 1, hour(t) + 1, t
end

function _checked(p, t)
    i = findfirst(isnan, p)
    isnothing(i) || throw(ArgumentError("OMNI has no `$i` at $t (no data or a gap over 3 hours)"))
    return p
end

function TsyganenkoModels.inputs(::OMNI, ::Type{<:T89}, t)
    d, _, h, t = _sample(t)
    kp10 = _checked((; kp = d.kp[h]), t).kp  # OMNI stores Kp × 10: 1-, 1, 1+ are 7, 10, 13
    return (; iopt = min(round(Int, kp10 / 10) + 1, 7))
end

function TsyganenkoModels.inputs(::OMNI, ::Type{<:T96}, t)
    d, _, h, t = _sample(t)
    pdyn, byimf, bzimf = d.sw_hourly[h]
    return _checked((; pdyn, dst = d.dst[h], byimf, bzimf), t)
end

function TsyganenkoModels.inputs(::OMNI, ::Type{<:T01}, t)
    d, i, _, t = _sample(t)
    g1, g2 = d.g[i]
    return _checked((; pdyn = d.pdyn[i], dst = d.symh[i], byimf = d.by[i], bzimf = d.bz[i], g1, g2), t)
end

function TsyganenkoModels.inputs(::OMNI, ::Type{<:TS04}, t)
    d, i, _, t = _sample(t)
    w1, w2, w3, w4, w5, w6 = d.w[i]
    return _checked((; pdyn = d.pdyn[i], dst = d.symh[i], byimf = d.by[i], bzimf = d.bz[i], w1, w2, w3, w4, w5, w6), t)
end

end
