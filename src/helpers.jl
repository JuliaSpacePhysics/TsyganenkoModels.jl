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
        sinpsas = sps / (1.0 + (r / rh)^3)^0.33333333
        cospsas = sqrt(1.0 - sinpsas^2)
        xss_new = x * cospsas - z * sinpsas
        zss_new = x * sinpsas + z * cospsas
        abs(xss_new - xss) + abs(zss_new - zss) < 1.0e-6 && break
        xss, zss = xss_new, zss_new
    end
    return xss, zss
end
