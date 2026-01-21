#
# Tsyganenko 96 Magnetospheric Magnetic Field Model (T96_01)
#
# Based on the Fortran implementation by N.A. Tsyganenko
# https://geo.phys.spbu.ru/~tsyganenko/models/t96/T96.for
#
# References:
# - Tsyganenko, N.A., and D.P. Stern, "Modeling the global magnetic field of the
#   large-scale Birkeland current systems", J.Geophys.Res., v.101, p.27187-27198, 1996.
# - Tsyganenko, N.A., "Modeling the Earth's magnetospheric magnetic field confined
#   within a realistic magnetopause", J.Geophys.Res., v.100, pp.5599-5612, 1995.
# - https://github.com/tsssss/geopack/blob/master/geopack/t96.py

include("t96_helpers.jl")
include("t96_consts.jl")

"""
    t96(x, y, z, ps, pdyn, dst, byimf, bzimf) -> (Bx, By, Bz)

Compute GSM components of the external magnetic field using the Tsyganenko 96 model.

# Parameters
- `x, y, z`: Position in GSM coordinates [Earth radii]
- `ps`: Geodipole tilt angle [radians]
- `pdyn`: Solar wind dynamic pressure [nPa]
- `dst`: Dst index [nT]
- `byimf`: IMF By component [nT]
- `bzimf`: IMF Bz component [nT]

# Returns
- `(Bx, By, Bz)`: Magnetic field components in GSM coordinates [nT]

# Model Description
Data-based model calibrated by solar wind pressure, Dst index, and IMF By/Bz components.
Includes realistic magnetopause, Region 1 and 2 Birkeland current systems, and IMF penetration.

Valid parameter ranges (caution needed outside these ranges):
- Pdyn: 0.5 to 10 nPa
- Dst: -100 to +20 nT
- ByIMF and BzIMF: -10 to +10 nT

# References
- Tsyganenko & Stern, JGR, v.101, p.27187-27198, 1996
- Tsyganenko, JGR, v.100, p.5599-5612, 1995
"""
function t96(x, y, z, ps, pdyn, dst, byimf, bzimf)
    sps = sin(ps)
    depr = 0.8 * dst - 13.0 * sqrt(pdyn)

    bt = sqrt(byimf^2 + bzimf^2)
    if byimf == 0.0 && bzimf == 0.0
        theta = 0.0
    else
        theta = atan(byimf, bzimf)
        if theta < 0.0
            theta += 2 * π
        end
    end
    st, ct = sincos(theta)
    eps = 718.5 * sqrt(pdyn) * bt * sin(theta / 2)

    facteps = eps / EPS10 - 1.0
    factpd = sqrt(pdyn / PDYN0) - 1.0
    rcampl = -A_COEF[1] * depr

    tampl2 = A_COEF[2] + A_COEF[3] * factpd + A_COEF[4] * facteps
    tampl3 = A_COEF[5] + A_COEF[6] * factpd
    b1ampl = A_COEF[7] + A_COEF[8] * facteps
    b2ampl = 20.0 * b1ampl
    reconn = A_COEF[9]

    xappa = (pdyn / PDYN0)^0.14
    xappa3 = xappa^3
    ys = y * ct - z * st
    zs = z * ct + y * st

    factimf = exp(x / DELIMFX - (ys / DELIMFY)^2)
    oimfx = 0.0
    oimfy = reconn * byimf * factimf
    oimfz = reconn * bzimf * factimf
    oimf = (oimfx, oimfy, oimfz)

    rimfampl = reconn * bt

    xx = x * xappa
    yy = y * xappa
    zz = z * xappa

    x0 = X00 / xappa
    am = AM0 / xappa
    rho2 = y^2 + z^2
    asq = am^2
    xmxm = am + x - x0
    if xmxm < 0
        xmxm = 0.0
    end
    axx0 = xmxm^2
    aro = asq + rho2
    sigma = sqrt((aro + axx0 + sqrt((aro + axx0)^2 - 4 * asq * axx0)) / (2 * asq))

    return if sigma < (S0 + DSIG)
        cf = dipshld(ps, xx, yy, zz)
        brc, bt2, bt3 = tailrc96(sps, xx, yy, zz)
        r1 = birk1tot_02(ps, xx, yy, zz)
        r2 = birk2tot_02(ps, xx, yy, zz)
        rimfx, rimfys, rimfzs = intercon(xx, ys * xappa, zs * xappa)
        rimfy = rimfys * ct + rimfzs * st
        rimfz = rimfzs * ct - rimfys * st
        f = @. (rimfx, rimfy, rimfz) * rimfampl + cf * xappa3 + b1ampl * r1 + b2ampl * r2 + rcampl * brc + tampl2 * bt2 + tampl3 * bt3

        return if sigma < (S0 - DSIG)
            f
        else
            fint = 0.5 * (1 - (sigma - S0) / DSIG)
            fext = 0.5 * (1 + (sigma - S0) / DSIG)
            q = dipole(ps, x, y, z)
            @. (f + q) * fint + oimf * fext - q
        end
    else
        oimf .- dipole(ps, x, y, z)
    end
end

# Allow vector input for convenience
@inline function t96(r, args...)
    @assert length(r) == 3
    return t96(r[1], r[2], r[3], args...)
end
