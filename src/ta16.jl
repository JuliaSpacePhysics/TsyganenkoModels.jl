# Tsyganenko-Andreeva 2016 RBF model (TA16), from TA16_RBF.f (corrected version of 11/23/2021) and its
# coefficient file TA16_RBF.par, https://geo.phys.spbu.ru/~tsyganenko/empirical-models/magnetic_field/ta16
# (unreachable 2026-10; source archived at
# https://web.archive.org/web/20230922185916/https://geo.phys.spbu.ru/~tsyganenko/models/ta16_rbf/TA16_RBF.f,
# coefficients mirrored in https://github.com/drsteve/LANLGeoMag and https://github.com/sastanish/magnetosphere-T-models).

"""
    TA16(; pdyn, symhc, xind, byimf, vgse = nothing)

Tsyganenko-Andreeva 2016 RBF model: the field of 1296 radial basis functions and their mirror images in the
SM equatorial plane, warped with dipole tilt.

The model is fitted in the aberrated GSW frame, `GSW(vgse)`; `vgse = nothing` takes GSM.

# Parameters
- `pdyn`: Solar wind dynamic pressure [nPa]
- `symhc`: Corrected SYM-H [nT], the 30-min centered sliding average of `0.8 SYM-H - 13 √pdyn`
- `xind`: Newell et al. (2007) coupling index as normalized by Tsyganenko & Andreeva (2015), averaged over the
  preceding 30 min; typically 0 (quiet) to 2 (strongly disturbed)
- `byimf`: IMF By [nT] in the model frame, averaged over the preceding 30 min
- `vgse`: solar wind velocity relative to Earth (GSE) [km/s], see `GSW`

Valid only sunward of X = -15 Re.

# References
- Tsyganenko, N. A., and Andreeva, V. A. (2016), An empirical RBF model of the magnetosphere parameterized by
  interplanetary and ground-based drivers, J. Geophys. Res. Space Physics, https://doi.org/10.1002/2016JA023217.
- Andreeva, V. A., and Tsyganenko, N. A. (2016), Reconstructing the magnetosphere from data using radial basis
  functions, J. Geophys. Res. Space Physics, https://doi.org/10.1002/2015JA022242.
"""
struct TA16{T, F} <: TsyganenkoModel
    pdyn::T
    symhc::T
    xind::T
    byimf::T
    frame::F
end

TA16(; pdyn, symhc, xind, byimf, vgse = nothing) =
    TA16(promote(pdyn, symhc, xind, byimf)..., isnothing(vgse) ? GSM() : GSW(vgse))

getcsys(m::TA16) = (m.frame, Cartesian3())

const TA16_N = 1296
const TA16_D² = 16.0  # squared RBF smoothing scale [Re²]

# Coefficients of center i: the toroidal-cos, toroidal-sin, poloidal-cos, poloidal-sin amplitudes as
# 4×4 (amplitude × (1, fpd, symh, xind)), then the IMF By toroidal and poloidal amplitudes.
const TA16_COEF = let path = joinpath(@__DIR__, "..", "data", "TA16_RBF.par")
    include_dependency(path)
    a = parse.(Float64, Iterators.take(eachline(path), 18TA16_N))
    [ntuple(k -> a[i + TA16_N * (k - 1)], 18) for i in 1:TA16_N]
end

# Lin et al. (2010) magnetopause at Pdyn = 0.5 nPa, IMF Bz = 0 and zero tilt, at angle θ from the X axis: RBF centers beyond it are not in the grid.
function ta16_magnetopause(θ)
    s, c = sincos(θ)
    return 12.544 * 0.5^-0.194 * (sqrt((1 + c) / 2) + 0.0571 * 2s * c * (1 - exp(-θ)))^-0.999
end

# RBF centers on a Kurihara grid in the northern SM hemisphere, ordered as the coefficients: spheres growing
# from r = 3.3 Re, then rings of colatitude θ from the pole (none on the equator), then longitude φ. Per ring:
# `warp`, the tilt hinging factor of the sphere (Tsyganenko & Andreeva 2015, RH = 8, α = 3), and the range of
# its centers in `TA16_SINCOSφ`.
const TA16_RINGS, TA16_SINCOSφ = let nlat = 9, rh = 8.0
    rings = @NamedTuple{r::Float64, θ::Float64, warp::Float64, idx::UnitRange{Int}}[]
    sincosφ = NTuple{2, Float64}[]
    r = 3.3
    while true
        warp = rh / r * (1 - cbrt(1 + (r / rh)^3))
        for i in 1:nlat
            θ = (i - 1) * π / 2 / (nlat - 0.5)
            i0 = length(sincosφ) + 1
            nlon = max(4(i - 1), 1)
            for k in 1:nlon
                sφ, cφ = sincospi(2(k - 1) / nlon)
                x, ρ = r * sin(θ) * cφ, r * hypot(sin(θ) * sφ, cos(θ))
                r ≤ ta16_magnetopause(atan(ρ, x)) && push!(sincosφ, (sφ, cφ))
            end
            push!(rings, (; r, θ, warp, idx = i0:length(sincosφ)))
        end
        r *= (nlat - 0.5 + π / 4) / (nlat - 0.5 - π / 4)
        r > 16 && break
    end
    @assert length(sincosφ) == TA16_N
    rings, sincosφ
end

function evalmodel(m::TA16, x, y, z, ps::Real)
    (; pdyn, symhc, xind, byimf) = m
    sps, cps = sincos(ps)
    tps = sps / cps
    # RBF expansions are in SM coordinates
    xs, ys, zs = x * cps - z * sps, y, z * cps + x * sps
    fpd = sqrt(pdyn / 2) - 1
    symh = symhc / 50
    T = float(promote_type(typeof(xs), typeof(fpd)))
    r = SVector(xs, ys, zs)
    acc = (zero(SVector{3, T}), zero(SVector{3, T}), zero(T))
    for ring in TA16_RINGS
        # Centers shift in colatitude with tilt; mirror images (θ → π - θ) shift the same way.
        dθ = -asin(ring.warp * tps) * sin(ring.θ)
        sn, cn = ring.r .* sincos(ring.θ + dθ)
        ss, cs = ring.r .* sincos(ring.θ - dθ)
        for i in ring.idx
            sφ, cφ = TA16_SINCOSφ[i]
            a = TA16_COEF[i]
            amp(j) = a[j] + a[j + 4] * fpd + a[j + 8] * symh + a[j + 12] * xind
            tc, ts, pc, psn = cps * amp(1), sps * amp(2), cps * amp(3), sps * amp(4)
            tb, pb = a[17] * byimf, a[18] * byimf
            acc = ta16_rbf(acc, r, SVector(sn * cφ, sn * sφ, cn), tc + ts + tb, pc + psn + pb)
            acc = ta16_rbf(acc, r, SVector(ss * cφ, ss * sφ, -cs), tc - ts - tb, psn + pb - pc)
        end
    end
    v, u, s = acc
    bx, by, bz = r × v + s * r - u
    return SVector(bx * cps + bz * sps, by, bz * cps - bx * sps)
end

# Add the toroidal field ∇C × r and poloidal field ∇ × (∇C × r), weighted by `wt` and `wp`, of the RBF
# C = √(|r - c|² + D²) centered at c to the accumulators (v, u, s) of B = r × v + s r - u:
# ∇C × r = r × c / C, and ∇ × (∇C × r) = [(1 - (d⋅r + D²)/C²) r - (2 - d⋅r/C²) c] / C with d = r - c.
@inline function ta16_rbf((v, u, s), r, c, wt, wp)
    d = r - c
    q = 1 / (d ⋅ d + TA16_D²)
    ic = sqrt(q)
    dr = d ⋅ r
    p = wp * ic
    return (v + wt * ic * c, u + p * (2 - dr * q) * c, s + p * (1 - (dr + TA16_D²) * q))
end
