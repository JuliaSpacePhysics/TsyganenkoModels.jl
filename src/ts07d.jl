# Tsyganenko-Sitnov TS07D model, from TS07D_JULY_2017 (https://github.com/grantkstephens/ARGUS, src/ts07, whose June
# 2023 restructuring is the same model; mirrored with modern-compiler fixes in https://github.com/scivision/ts07d)
# and its static tail shielding coefficients (TAIL_PAR.zip of https://doi.org/10.5281/zenodo.14733741,
# data/TS07D_TAIL_PAR).

"""
    TS07D(; a, pdyn)
    TS07D(file)

Tsyganenko-Sitnov TS07D model with its 101 fitted coefficients and parameters `a` and solar wind dynamic pressure
`pdyn` [nPa], or both from a coefficient `file` of the TS07D archive (`YYYY_DDD_HH_MM.par`, one per 5 min for
1995–2018, at https://doi.org/10.5281/zenodo.14733741).

The equatorial currents are expanded in 5 radial × (1 + 2×4 azimuthal) Bessel-function modes, with dipole shielding
and region-1 and -2 Birkeland currents; there is no magnetopause, so the field is only meaningful inside it.

# References
- Tsyganenko, N. A., and Sitnov, M. I. (2007), Magnetospheric configurations from a high-resolution data-based
  magnetic field model, J. Geophys. Res., 112, A06225, https://doi.org/10.1029/2007JA012260.
- Sitnov, M. I., et al. (2008), Dynamical data-based modeling of the storm-time geomagnetic field with enhanced
  spatial resolution, J. Geophys. Res., 113, A07218, https://doi.org/10.1029/2007JA013003.
"""
struct TS07D{T} <: TsyganenkoModel
    a::NTuple{101, T}
    pdyn::T
end

function TS07D(; a, pdyn)
    T = promote_type(eltype(a), typeof(pdyn))
    return TS07D{T}(NTuple{101, T}(a), pdyn)
end

# 101 lines of a value and a comment, then `Q=`, `B_rms=`, `M=`, `N=`, `Pdyn=`, `tilt=` lines
function TS07D(file::AbstractString)
    lines = readlines(file)
    a = ntuple(i -> parse(Float64, first(split(lines[i]))), 101)
    pdyn = parse(Float64, split(only(filter(startswith("Pdyn="), strip.(lines))), '=')[2])
    return TS07D(; a, pdyn)
end

evalmodel(m::TS07D, x, y, z, ps::Real) = GSM(ts07d_extern(m.a, m.pdyn, ps, x, y, z))

# Shielding of the tail modes: per mode, 75 coefficients c[n + 5m] (n = 1:5, m = 0:14) and wave numbers |c[75 + n]|;
# symmetric modes k, odd and even modes (k, l)
const TS07D_TSS, TS07D_TSO, TS07D_TSE = let dir = joinpath(@__DIR__, "..", "data", "TS07D_TAIL_PAR")
    function load(f)
        path = joinpath(dir, f)
        include_dependency(path)
        lines = readlines(path)
        return ntuple(i -> parse(Float64, lines[i]), 80)
    end
    (
        ntuple(k -> load("tailamebhr$k.par"), 5),
        ntuple(k -> ntuple(l -> load("tailamhr_o_$k$l.par"), 4), 5),
        ntuple(k -> ntuple(l -> load("tailamhr_e_$k$l.par"), 4), 5),
    )
end

function ts07d_extern(a, pdyn, ps, x, y, z)
    xappa = (pdyn / 2)^0.155
    xx, yy, zz = x * xappa, y * xappa, z * xappa
    d, rh0, g, xkappa1, xkappa2, tw = a[96], a[97], a[98], a[99], a[100], a[101]

    pf = sqrt(pdyn / 2) - 1
    ws = ntuple(k -> a[1 + k] + a[46 + k] * pf, 5)
    wo = ntuple(k -> ntuple(l -> a[6 + 4(k - 1) + l] + a[51 + 4(k - 1) + l] * pf, 4), 5)
    we = ntuple(k -> ntuple(l -> a[26 + 4(k - 1) + l] + a[71 + 4(k - 1) + l] * pf, 4), 5)
    # the deformations are linear in the field, so the weighted sum of the modes is deformed once
    (bt,) = deformed((x, y, z) -> (ts07d_unwarped(x, y, z, d, ws, wo, we),), ps, xx, yy, zz, rh0, g, tw)

    bcf = dipole_shield(xx, yy, zz, ps) .* xappa^3
    br11, br12, br21a, _ = birk_tot(ps, xx, yy, zz, xkappa1, xkappa2)
    br21s = birk_1n2(2, 1, ps, xx, yy, zz, xkappa2, twoconss) .+
        _harmonic_shield(SH21S, ps, xkappa2 - 1.0, xx, yy, zz, -1.0, true)
    return a[1] .* bcf .+ bt .+ a[92] .* br11 .+ a[93] .* br12 .+ a[94] .* br21a .+ a[95] .* br21s
end

# Weighted sum of the shielded tail modes (sheet thickness `d`) before warping and bending
function ts07d_unwarped(x, y, z, d, ws, wo, we)
    rho = sqrt(x^2 + y^2)
    scφ = ntuple(m -> sincos((m - 1) * atan(y, x)), Val(15))  # sin, cos of mφ for m = 0:14
    sφ, cφ = scφ[2]
    zd = sqrt(z^2 + d^2)
    bx = by = bz = zero(rho)
    for k in 1:5
        rkm = k / 20
        rkmr = rkm * rho
        J = besselj_upto(Val(5), rkmr)
        rex = exp(rkm * zd)
        t = ws[k] * rkm * z * J[2] / zd / rex
        bx += t * cφ
        by += t * sφ
        bz += ws[k] * rkm * J[1] / rex
        for m in 1:4
            Jm = J[m + 1]
            Jx, dJm = besselj_over_x(J, m, rkmr)
            sm, cm = scφ[m + 1]
            o, e = wo[k][m], we[k][m]
            bρ = m * z * dJm / zd / rex * (o * cm - e * sm)
            bϕ = -m^2 * z * Jx / zd / rex * (o * sm + e * cm)
            bx += bρ * cφ - bϕ * sφ
            by += bρ * sφ + bϕ * cφ
            bz += m * Jm / rex * (e * sm - o * cm)
        end
        bx, by, bz = (bx, by, bz) .+ ws[k] .* bessel_shield(TS07D_TSS[k], false, scφ, z, rho)
        for l in 1:4
            bx, by, bz = (bx, by, bz) .+
                wo[k][l] .* bessel_shield(TS07D_TSO[k][l], false, scφ, z, rho) .+
                we[k][l] .* bessel_shield(TS07D_TSE[k][l], true, scφ, z, rho)
        end
    end
    return bx, by, bz
end

# -∇ of Σₙ Σₘ c[n + 5m] T(mφ) Jₘ(kₙρ) sinh(kₙz), kₙ = |c[75 + n]|, with T = sin if `odd`, else cos
function bessel_shield(c, odd, scφ, z, rho)
    sφ, cφ = scφ[2]
    hx = hy = hz = zero(rho)
    for n in 1:5
        k = abs(c[75 + n])
        kρ = k * rho
        shz, chz = sinh(z * k), cosh(z * k)
        J = besselj_upto(Val(14), kρ)
        for m in 0:14
            sm, cm = scφ[m + 1]
            T, dT = odd ? (sm, m * cm) : (cm, -m * sm)  # T(mφ) and its φ derivative
            Jx, dJm = besselj_over_x(J, m, kρ)
            cnm = c[n + 5m]
            # ∇ρ = (cos φ, sin φ), ρ∇φ = (-sin φ, cos φ)
            hx -= cnm * k * shz * (T * dJm * cφ - dT * Jx * sφ)
            hy -= cnm * k * shz * (T * dJm * sφ + dT * Jx * cφ)
            hz -= cnm * k * T * chz * J[m + 1]
        end
    end
    return hx, hy, hz
end

# J₀(x), …, J_N(x) for x ≥ 0 by Miller's downward recurrence, normalized by J₀ + 2(J₂ + J₄ + …) = 1, as TS07D's `bessJJ`
function besselj_upto(::Val{N}, x) where {N}
    # Leading series term, exact to double precision; the recurrence's growth (~2j/x per step) would overflow
    x < 1.0e-8 && return ntuple(i -> (x / 2)^(i - 1) / factorial(i - 1), Val(N + 1))
    J = ntuple(_ -> zero(x), Val(N + 1))
    tox = 2 / x
    evnsum, iseven = zero(x), false
    bjp, bj = zero(x), one(x)
    for j in (2 * ((N + isqrt(40N)) ÷ 2)):-1:1
        bjp, bj = bj, j * tox * bj - bjp
        if abs(bj) > 1.0e10
            bj, bjp, evnsum = bj * 1.0e-10, bjp * 1.0e-10, evnsum * 1.0e-10
            J = _scaled_above(J, j + 1, 1.0e-10)
        end
        iseven && (evnsum += bj)
        iseven = !iseven
        j <= N && (J = _setindex(J, bjp, j + 1))
    end
    return _setindex(J, bj, 1) ./ (2evnsum - bj)
end

# Jₘ(x)/x and Jₘ'(x) from `J` = (J₀, J₁, …), finite on the axis x = 0. Jₘ/x is only needed for m ≥ 1, as it
# multiplies the φ derivative of the mode.
@inline function besselj_over_x(J, m, x)
    m == 0 && return zero(x), -J[2]
    Jx = x == 0 ? (m == 1 ? one(x) / 2 : zero(x)) : J[m + 1] / x
    return Jx, J[m] - m * Jx
end

# Unrolled tuple updates: `Base.setindex` on long tuples with a runtime index does not infer
_setindex(t::NTuple{N}, v, k) where {N} = ntuple(i -> i == k ? v : t[i], Val(N))
_scaled_above(t::NTuple{N}, k, s) where {N} = ntuple(i -> i > k ? t[i] * s : t[i], Val(N))

const SH21S = [
    -67763516.61, -49565522.84, 10123356.08, 51805446.1, -51607711.68, 164360662.1,
    -4662006.024, -191297217.6, -7204547.103, 30372354.93, -750371.9365, -36564457.17,
    61114395.65, 45702536.5, -9228894.939, -47893708.68, 47290934.33, -149155112.0,
    4226520.638, 173588334.5, 7998505.443, -33150962.72, 832493.2094, 39892545.84,
    -11303915.16, -8901327.398, 1751557.11, 9382865.82, -9054707.868, 27918664.5,
    -788741.7146, -32481294.42, -2264443.753, 9022346.503, -233526.0185, -10856269.53,
    -244450.885, 1908295.272, 185445.1967, -1074202.863, 41827.75224, -241553.7626,
    -20199.1258, 123235.6084, 199501.4614, -1936498.464, -178857.4074, 1044724.507,
    121044.9917, -946479.9247, -91808.28803, 532742.7569, -20742.28628, 120633.2193,
    10018.49534, -61599.11035, -98709.58977, 959095.177, 88500.43489, -517471.5287,
    -81.56122911, 816.2472344, 55.3071171, -454.5368824, 25.7469381, -202.500735,
    -7.369350794, 104.9429812, 58.14049362, -685.5919355, -51.71345683, 374.0125033,
    247.9296982, 159.2471769, 102.3151816, 15.81062488, 34.99767599, 133.0832773,
    219.6475201, 107.9582783, 10.00264684, 7.718306072, 25.22866153, 5.013583103,
    0.08407754233, -9.613356793e-06,
]
