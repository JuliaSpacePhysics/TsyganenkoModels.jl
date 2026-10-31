module TsyganenkoModelsCDAWebExt

using CDAWeb: CDAWeb, getdata, times
using Dates
using GeoCotrans: GSM, GSW, SVector, rotation
using TsyganenkoModels: TsyganenkoModels, OMNI, T89, T96, T01, TS04, TA16, fill_gaps!, masked_mean, t01_g, ta16_inputs, ts05_w, TS05_R

const STEP = Minute(5)
const PER_DAY = 288
# 5-min samples fetched before the day, for W memory to decay below e⁻¹⁰ at the slowest rate (W3, ~13 days)
const NBEFORE = ceil(Int, 10 / (minimum(TS05_R) * 5 / 60))
# Samples kept on each side of the day: the widest window, T01's 12 samples ending at `t`
const PAD = 12
# A bounding observation beyond the longest interpolated run (35 missing samples).
const GAP_PAD = 36

# Gap-filled 5-min observations shared by the models, including PAD samples on each side of the day
struct OMNIDay
    pdyn::Vector{Float64}
    v::Vector{Float64}
    # IMF in GSM; OMNI has only BX_GSE = GSM Bx
    bx::Vector{Float64}
    by::Vector{Float64}
    bz::Vector{Float64}
    symh::Vector{Float64}
end

# On the OMNIDay grid; IMF in the GSW frame of each sample's `vgse`
struct TA16Obs
    by::Vector{Float64}
    bz::Vector{Float64}
    vgse::Vector{NTuple{3, Float64}}
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

function _fetch_5min(names, t0, n)
    hro = CDAWeb.Dataset("OMNI_HRO_5MIN")
    return map(names) do name
        fill_gaps!(_series(hro, name, t0, n, STEP))
    end
end

# Gap-filled series of `names` on the OMNIDay grid
function _fetch_5min(names, day::Date)
    margin = PAD + GAP_PAD
    xs = _fetch_5min(names, DateTime(day) - margin * STEP, PER_DAY + 2margin)
    return map(x -> x[(GAP_PAD + 1):(end - GAP_PAD)], xs)
end

function OMNIDay(day::Date)
    den, v, bx, by, bz, symh = _fetch_5min(("proton_density", "flow_speed", "BX_GSE", "BY_GSM", "BZ_GSM", "SYM_H"), day)
    pdyn = @. 1.937e-6 * den * v^2  # 1.16 mₚ in nPa per (cm⁻³ (km/s)²), Tsyganenko's value
    return OMNIDay(pdyn, v, bx, by, bz, symh)
end

function TA16Obs(day::Date)
    d = _omni_day(day)
    vx, vy, vz = _fetch_5min(("Vx", "Vy", "Vz"), day)
    vgse = tuple.(vx, vy .+ 29.78, vz)
    t0 = DateTime(day) - PAD * STEP
    b = map(eachindex(vgse)) do k
        rotation(GSM() => GSW(vgse[k]), t0 + (k - 1) * STEP) * SVector(d.bx[k], d.by[k], d.bz[k])
    end
    return TA16Obs(getindex.(b, 2), getindex.(b, 3), vgse)
end

function _ts04_w(day)
    t0 = DateTime(day) - NBEFORE * STEP
    n = NBEFORE + PER_DAY + GAP_PAD
    den, v, bz = _fetch_5min(("proton_density", "flow_speed", "BZ_GSM"), t0, n)
    return ts05_w(den, v, bz)[NBEFORE .+ (1:PER_DAY)]
end

const LOCK = ReentrantLock()
const OMNI_DAYS = Dict{Date, OMNIDay}()
const HOURLY = Dict{Tuple{Date, String}, Vector{Float64}}()
const TS04_W = Dict{Date, Vector{NTuple{6, Float64}}}()
const TA16_OBS = Dict{Date, TA16Obs}()

_cached(f, cache, key) = @lock LOCK get!(f, cache, key)
_omni_day(day) = _cached(() -> OMNIDay(day), OMNI_DAYS, day)
_hourly(name, day) = _cached(HOURLY, (day, name)) do
    _series(CDAWeb.Dataset("OMNI2_H0_MRG1HR"), name, DateTime(day), 24, Hour(1))
end

function _checked(p, t)
    i = findfirst(x -> any(isnan, x), p)
    isnothing(i) || throw(ArgumentError("OMNI has no `$i` at $t (no data or a gap over 3 hours)"))
    return p
end

function TsyganenkoModels.inputs(::OMNI, M::Type{<:Union{T89, T96, T01, TS04, TA16}}, t)
    t = DateTime(t)
    day = Date(t)
    return _input(M, day, (t - DateTime(day)) ÷ STEP + 1, t)
end

# Sample `j` of `day` is `PAD + j` in its `OMNIDay`
function _input(::Type{<:T89}, day, _, t)
    kp10 = _checked((; kp = _hourly("KP", day)[hour(t) + 1]), t).kp  # OMNI stores Kp × 10: 1-, 1, 1+ are 7, 10, 13
    return (; iopt = min(round(Int, kp10 / 10) + 1, 7))
end

function _input(::Type{<:T96}, day, _, t)
    d, h = _omni_day(day), hour(t) + 1
    pdyn, byimf, bzimf = map(x -> masked_mean(view(x, (PAD + 12h - 11):(PAD + 12h))), (d.pdyn, d.by, d.bz))
    return _checked((; pdyn, dst = _hourly("DST", day)[h], byimf, bzimf), t)
end

function _input(::Type{<:T01}, day, j, t)
    d, i = _omni_day(day), PAD + j
    g1, g2 = t01_g(d.v, d.by, d.bz, i)
    return _checked((; pdyn = d.pdyn[i], dst = d.symh[i], byimf = d.by[i], bzimf = d.bz[i], g1, g2), t)
end

function _input(::Type{<:TS04}, day, j, t)
    d, i = _omni_day(day), PAD + j
    w1, w2, w3, w4, w5, w6 = _cached(() -> _ts04_w(day), TS04_W, day)[j]
    return _checked((; pdyn = d.pdyn[i], dst = d.symh[i], byimf = d.by[i], bzimf = d.bz[i], w1, w2, w3, w4, w5, w6), t)
end

function _input(::Type{<:TA16}, day, j, t)
    d, o, i = _omni_day(day), _cached(() -> TA16Obs(day), TA16_OBS, day), PAD + j
    return _checked((; pdyn = d.pdyn[i], ta16_inputs(d.pdyn, d.symh, d.v, o.by, o.bz, i)..., vgse = o.vgse[i]), t)
end

end
