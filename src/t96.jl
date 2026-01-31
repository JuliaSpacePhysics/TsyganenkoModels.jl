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
    x0 = X00 / xappa
    am = AM0 / xappa
    rho2 = y^2 + z^2
    sigma = _sigma(x, x0, am, rho2)

    s0 = 1.08
    dsig = 0.005
    out = _switch(sigma, s0, dsig, ps, x, y, z, oimf; q0 = 30574) do
        xx = x * xappa; yy = y * xappa; zz = z * xappa
        cf = dipshld(ps, xx, yy, zz)
        brc, bt2, bt3 = tailrc96(sps, xx, yy, zz)
        r1 = birk1tot_02(ps, xx, yy, zz)
        r2 = birk2tot_02(ps, xx, yy, zz)
        rimfx, rimfys, rimfzs = intercon(xx, ys * xappa, zs * xappa)
        rimfy = rimfys * ct + rimfzs * st
        rimfz = rimfzs * ct - rimfys * st
        rimf = (rimfx, rimfy, rimfz)
        @. rimf * rimfampl + cf * xappa3 + b1ampl * r1 + b2ampl * r2 + rcampl * brc + tampl2 * bt2 + tampl3 * bt3
    end
    return GSM(out)
end

function _sigma(x, x0, am, rho2)
    xmxm = max(am + x - x0, 0.0)
    asq = am^2
    aro = asq + rho2
    axx0 = xmxm^2
    return sqrt((aro + axx0 + sqrt((aro + axx0)^2 - 4 * asq * axx0)) / (2 * asq))
end

evalmodel(m::T96, x, y, z, ps) = t96(x, y, z, ps, m.pdyn, m.dst, m.byimf, m.bzimf)
