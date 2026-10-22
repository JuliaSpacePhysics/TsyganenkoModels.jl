# Shue et al. magnetopause parameters, shared by T01 and TS04.
const MP_A  = 34.586
const MP_S0 = 1.196
const MP_X0 = 3.4397

@inline function _switch(f, sigma, s0, dsig, ps, x, y, z, oimf; q0)
    return if sigma < (s0 + dsig)
        b = f()
        if sigma < (s0 - dsig)
            b
        else # Interpolation region near magnetopause
            fint = 0.5 * (1 - (sigma - s0) / dsig)
            fext = 1.0 - fint
            q = dipole(ps, x, y, z; q0)
            @. (b + q) * fint + oimf * fext - q
        end
    else  # Outside magnetosphere
        oimf .- dipole(ps, x, y, z; q0)
    end
end

# Fixed-point iteration for tail-warped coordinates (xss, zss).
# Solves the Tsyganenko hinging equation by Picard iteration; converges in <20 steps.
@inline function _warped_xz(x, y, z, sps, rh0; rh2 = -5.2)
    r = sqrt(x^2 + y^2 + z^2)
    xss, zss = x, z
    for _ in 1:20
        rh = rh0 + rh2 * (zss / r)^2
        sinpsas = sps / cbrt(1.0 + (r / rh)^3)
        cospsas = sqrt(1.0 - sinpsas^2)
        xss_new = x * cospsas - z * sinpsas
        zss_new = x * sinpsas + z * cospsas
        abs(xss_new - xss) + abs(zss_new - zss) < 1.0e-6 && break
        xss, zss = xss_new, zss_new
    end
    return xss, zss
end

"""
    cartharm(c, p, r, x, y, z, parallel) -> (Bx, By, Bz)

Evaluate the Cartesian-harmonic field `B = Σ c(i,k) ∇Uᵢₖ`, where
`Uᵢₖ = -exp(x√(1/pᵢ² + 1/rₖ²)) cos(y/pᵢ) f(z/rₖ)`.
`p` and `r` are spatial scales along y and z; `parallel` selects `f = cos`
(Bx even, Bz odd in z), otherwise `f = sin` (Bx odd, Bz even).
Each potential satisfies Laplace's equation. Fitted coefficients provide magnetopause
shielding or IMF penetration (Tsyganenko 1995, JGR 100, 5599).
"""
@inline function cartharm(c, p, r, x, y, z, parallel)
    syp = map(pᵢ -> sincos(y / pᵢ), p)
    szr = map(rₖ -> sincos(z / rₖ), r)
    bx = by = bz = 0.0
    for i in eachindex(p), k in eachindex(r)
        sy, cy = syp[i]
        sz, cz = szr[k]
        fz, dfz = parallel ? (cz, -sz) : (sz, cz)
        sq = sqrt(1 / p[i]^2 + 1 / r[k]^2)
        e = c(i, k) * exp(x * sq)
        bx -= sq * e * cy * fz
        by += e / p[i] * sy * fz
        bz -= e / r[k] * cy * dfz
    end
    return bx, by, bz
end

# Field `f(x, y, z)` of a source tilted by `θ` about the y axis, as for the T01/TS04 shields.
@inline function tilted(f, θ, x, y, z)
    st, ct = sincos(θ)
    hx, hy, hz = f(x * ct - z * st, y, x * st + z * ct)
    return hx * ct + hz * st, hy, hz * ct - hx * st
end
