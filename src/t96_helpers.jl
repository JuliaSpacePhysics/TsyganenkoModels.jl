using Bumper
using SpecialFunctions

function dipshld(ps, x, y, z)
    sps, cps = sincos(ps)
    return cps .* cylharm(A1_DIPSHLD, x, y, z) .+ sps .* cylhar1(A2_DIPSHLD, x, y, z)
end

function cylharm(a, x, y, z)
    rho = sqrt(y^2 + z^2)
    if rho < 1.0e-8
        sinfi, cosfi = 1.0, 0.0
        rho = 1.0e-8
    else
        sinfi = z / rho
        cosfi = y / rho
    end
    sinfi2 = sinfi^2
    si2co2 = sinfi2 - cosfi^2

    bx = 0.0
    by = 0.0
    bz = 0.0
    for i in 1:3
        dzeta = rho / a[i + 6]
        xksi = x / a[i + 6]
        xj0 = besselj0(dzeta)
        xj1 = besselj1(dzeta)
        xexp = exp(xksi)
        bx -= a[i] * xj1 * xexp * sinfi
        by += a[i] * (2 * xj1 / dzeta - xj0) * xexp * sinfi * cosfi
        bz += a[i] * (xj1 / dzeta * si2co2 - xj0 * sinfi2) * xexp
    end
    for i in 4:6
        dzeta = rho / a[i + 6]
        xksi = x / a[i + 6]
        xj0 = besselj0(dzeta)
        xj1 = besselj1(dzeta)
        xexp = exp(xksi)
        brho = (xksi * xj0 - (dzeta^2 + xksi - 1) * xj1 / dzeta) * xexp * sinfi
        bphi = (xj0 + xj1 / dzeta * (xksi - 1)) * xexp * cosfi
        bx += a[i] * (dzeta * xj0 + xksi * xj1) * xexp * sinfi
        by += a[i] * (brho * cosfi - bphi * sinfi)
        bz += a[i] * (brho * sinfi + bphi * cosfi)
    end
    return bx, by, bz
end

function cylhar1(a, x, y, z)
    rho = sqrt(y^2 + z^2)
    if rho < 1.0e-8
        sinfi, cosfi = 1.0, 0.0
        rho = 1.0e-8
    else
        sinfi = z / rho
        cosfi = y / rho
    end

    bx = 0.0
    by = 0.0
    bz = 0.0
    for i in 1:3
        dzeta = rho / a[i + 6]
        xksi = x / a[i + 6]
        xj0 = besselj0(dzeta)
        xj1 = besselj1(dzeta)
        xexp = exp(xksi)
        brho = xj1 * xexp
        bx -= a[i] * xj0 * xexp
        by += a[i] * brho * cosfi
        bz += a[i] * brho * sinfi
    end
    for i in 4:6
        dzeta = rho / a[i + 6]
        xksi = x / a[i + 6]
        xj0 = besselj0(dzeta)
        xj1 = besselj1(dzeta)
        xexp = exp(xksi)
        brho = (dzeta * xj0 + xksi * xj1) * xexp
        bx += a[i] * (dzeta * xj1 - xj0 * (xksi + 1)) * xexp
        by += a[i] * brho * cosfi
        bz += a[i] * brho * sinfi
    end
    return bx, by, bz
end

# Calculates the potential interconnection field inside the magnetosphere, corresponding to
# DELTA_X = 20Re and DELTA_Y = 10Re (NB#3, p.90, 6/6/1996).
@views function intercon(x, y, z)
    a = INTERCON_A

    p = NTuple{3}(a[10:12])
    r = NTuple{3}(a[13:15])
    rp = 1.0 ./ p
    rr = 1.0 ./ r

    bx = 0.0
    by = 0.0
    bz = 0.0
    l = 1
    for i in 1:3
        sypi, cypi = sincos(y * rp[i])
        for k in 1:3
            szrk, czrk = sincos(z * rr[k])
            sqpr = sqrt(rp[i]^2 + rr[k]^2)
            epr = exp(x * sqpr)
            hx = -sqpr * epr * cypi * szrk
            hy = rp[i] * epr * sypi * szrk
            hz = -rr[k] * epr * cypi * czrk
            bx += a[l] * hx
            by += a[l] * hy
            bz += a[l] * hz
            l += 1
        end
    end
    return bx, by, bz
end

function shlcar3x3(a, x, y, z, sps)
    cps = sqrt(1 - sps^2)
    s3ps = 4 * cps^2 - 1

    hx = 0.0
    hy = 0.0
    hz = 0.0
    l = 1
    for m in 1:2
        for i in 1:3
            p = a[36 + i]
            q = a[42 + i]
            sypi, cypi = sincos(y / p)
            syqi, cyqi = sincos(y / q)
            for k in 1:3
                r = a[39 + k]
                s = a[45 + k]
                szrk, czrk = sincos(z / r)
                szsk, czsk = sincos(z / s)
                sqpr = sqrt(1 / p^2 + 1 / r^2)
                sqqs = sqrt(1 / q^2 + 1 / s^2)
                epr = exp(x * sqpr)
                eqs = exp(x * sqqs)
                if m == 1
                    dx = -sqpr * epr * cypi * szrk
                    dy = epr / p * sypi * szrk
                    dz = -epr / r * cypi * czrk
                    hx += a[l] * dx
                    hy += a[l] * dy
                    hz += a[l] * dz
                    l += 1

                    dx *= cps
                    dy *= cps
                    dz *= cps
                    hx += a[l] * dx
                    hy += a[l] * dy
                    hz += a[l] * dz
                    l += 1
                else
                    dx = -sps * sqqs * eqs * cyqi * czsk
                    dy = sps * eqs / q * syqi * czsk
                    dz = sps * eqs / s * cyqi * szsk
                    hx += a[l] * dx
                    hy += a[l] * dy
                    hz += a[l] * dz
                    l += 1

                    dx *= s3ps
                    dy *= s3ps
                    dz *= s3ps
                    hx += a[l] * dx
                    hy += a[l] * dy
                    hz += a[l] * dz
                    l += 1
                end
            end
        end
    end
    return hx, hy, hz
end

function tailrc96(sps, x, y, z)
    rh, dr, g, d0, deltady = 9.0, 4.0, 10.0, 2.0, 10.0
    dr2 = dr^2
    c11 = sqrt((1 + rh)^2 + dr2)
    c12 = sqrt((1 - rh)^2 + dr2)
    c1 = c11 - c12
    spsc1 = sps / c1
    rps = 0.5 * (c11 + c12) * sps

    r = sqrt(x^2 + y^2 + z^2)
    sq1 = sqrt((r + rh)^2 + dr2)
    sq2 = sqrt((r - rh)^2 + dr2)
    c = sq1 - sq2
    cs = (r + rh) / sq1 - (r - rh) / sq2
    spss = spsc1 / r * c
    cpss = sqrt(1 - spss^2)
    dpsrr = sps / (r * r) * (cs * r - c) / sqrt((r * c1)^2 - (c * sps)^2)

    wfac = y / (y^4 + 1.0e4)
    w = wfac * y^3
    ws = 4.0e4 * y * wfac^2
    warp = g * sps * w
    xs = x * cpss - z * spss
    zsww = z * cpss + x * spss
    zs = zsww + warp

    dxsx = cpss - x * zsww * dpsrr
    dxsy = -y * zsww * dpsrr
    dxsz = -spss - z * zsww * dpsrr
    dzsx = spss + x * xs * dpsrr
    dzsy = xs * y * dpsrr + g * sps * ws
    dzsz = cpss + xs * z * dpsrr

    d = d0 + deltady * (y / 20)^2
    dddy = deltady * y * 0.005

    dzetas = sqrt(zs^2 + d^2)
    ddzetadx = zs * dzsx / dzetas
    ddzetady = (zs * dzsy + d * dddy) / dzetas
    ddzetadz = zs * dzsz / dzetas

    warp_state = (;
        cpss, spss, dpsrr, rps, warp,
        d, xs, zs,
        dxsx, dxsy, dxsz, dzsx, dzsy, dzsz,
        dzetas, ddzetadx, ddzetady, ddzetadz,
        zsww,
    )

    brc = shlcar3x3(ARC_COEF, x, y, z, sps) .+ ringcurr96(x, y, z, warp_state)
    bt2 = shlcar3x3(ATAIL2_COEF, x, y, z, sps) .+ taildisk(x, y, z, warp_state)
    hx, hz = tail87(x, z, warp_state)
    bt3 = shlcar3x3(ATAIL3_COEF, x, y, z, sps) .+ (hx, 0, hz)
    return brc, bt2, bt3
end

function ringcurr96(x, y, z, warp)
    cpss = warp.cpss
    spss = warp.spss
    dpsrr = warp.dpsrr
    xs = warp.xs
    zsww = warp.zsww
    dxsx = warp.dxsx
    dxsy = warp.dxsy
    dxsz = warp.dxsz
    dzsx = warp.dzsx
    dzsz = warp.dzsz

    dzsy = xs * y * dpsrr
    xxd = x - XD_RC
    fdx = 0.5 * (1 + xxd / sqrt(xxd^2 + XLDX_RC^2))
    dddx = DELTADX_RC * 0.5 * XLDX_RC^2 / sqrt(xxd^2 + XLDX_RC^2)^3
    d = D0_RC + DELTADX_RC * fdx

    zs = zsww
    dzetas = sqrt(zs^2 + d^2)
    rhos = sqrt(xs^2 + y^2)
    ddzetadx = (zs * dzsx + d * dddx) / dzetas
    ddzetady = zs * dzsy / dzetas
    ddzetadz = zs * dzsz / dzetas

    if rhos < 1.0e-5
        drhosdx = 0.0
        drhosdy = sign(y)
        drhosdz = 0.0
    else
        drhosdx = xs * dxsx / rhos
        drhosdy = (xs * dxsy + y) / rhos
        drhosdz = xs * dxsz / rhos
    end

    bx = 0.0
    by = 0.0
    bz = 0.0
    for i in 1:2
        bi = BETA_RC[i]
        s1 = sqrt((dzetas + bi)^2 + (rhos + bi)^2)
        s2 = sqrt((dzetas + bi)^2 + (rhos - bi)^2)
        ds1ddz = (dzetas + bi) / s1
        ds2ddz = (dzetas + bi) / s2
        ds1drhos = (rhos + bi) / s1
        ds2drhos = (rhos - bi) / s2

        ds1dx = ds1ddz * ddzetadx + ds1drhos * drhosdx
        ds1dy = ds1ddz * ddzetady + ds1drhos * drhosdy
        ds1dz = ds1ddz * ddzetadz + ds1drhos * drhosdz

        ds2dx = ds2ddz * ddzetadx + ds2drhos * drhosdx
        ds2dy = ds2ddz * ddzetady + ds2drhos * drhosdy
        ds2dz = ds2ddz * ddzetadz + ds2drhos * drhosdz

        s1ts2 = s1 * s2
        s1ps2 = s1 + s2
        s1ps2sq = s1ps2^2
        fac1 = sqrt(s1ps2sq - (2 * bi)^2)
        as0 = fac1 / (s1ts2 * s1ps2sq)
        term1 = 1 / (s1ts2 * s1ps2 * fac1)
        fac2 = as0 / s1ps2sq
        dasds1 = term1 - fac2 / s1 * (s2^2 + s1 * (3 * s1 + 4 * s2))
        dasds2 = term1 - fac2 / s2 * (s1^2 + s2 * (3 * s2 + 4 * s1))

        dasdx = dasds1 * ds1dx + dasds2 * ds2dx
        dasdy = dasds1 * ds1dy + dasds2 * ds2dy
        dasdz = dasds1 * ds1dz + dasds2 * ds2dz

        bx += F_RC[i] * ((2 * as0 + y * dasdy) * spss - xs * dasdz + as0 * dpsrr * (y^2 * cpss + z * zs))
        by -= F_RC[i] * y * (as0 * dpsrr * xs + dasdz * cpss + dasdx * spss)
        bz += F_RC[i] * ((2 * as0 + y * dasdy) * cpss + xs * dasdx - as0 * dpsrr * (x * zs + y^2 * spss))
    end

    return bx, by, bz
end

function taildisk(x, y, z, warp)
    cpss = warp.cpss
    spss = warp.spss
    dpsrr = warp.dpsrr
    xs = warp.xs
    dxsx = warp.dxsx
    dxsy = warp.dxsy
    dxsz = warp.dxsz
    dzetas = warp.dzetas
    ddzetadx = warp.ddzetadx
    ddzetady = warp.ddzetady
    ddzetadz = warp.ddzetadz
    zsww = warp.zsww

    rhos = sqrt((xs - XSHIFT_TD)^2 + y^2)
    if rhos < 1.0e-5
        drhosdx = 0.0
        drhosdy = sign(y)
        drhosdz = 0.0
    else
        drhosdx = (xs - XSHIFT_TD) * dxsx / rhos
        drhosdy = ((xs - XSHIFT_TD) * dxsy + y) / rhos
        drhosdz = (xs - XSHIFT_TD) * dxsz / rhos
    end

    bx = 0.0
    by = 0.0
    bz = 0.0
    for i in 1:4
        bi = BETA_TD[i]
        s1 = sqrt((dzetas + bi)^2 + (rhos + bi)^2)
        s2 = sqrt((dzetas + bi)^2 + (rhos - bi)^2)
        ds1ddz = (dzetas + bi) / s1
        ds2ddz = (dzetas + bi) / s2
        ds1drhos = (rhos + bi) / s1
        ds2drhos = (rhos - bi) / s2

        ds1dx = ds1ddz * ddzetadx + ds1drhos * drhosdx
        ds1dy = ds1ddz * ddzetady + ds1drhos * drhosdy
        ds1dz = ds1ddz * ddzetadz + ds1drhos * drhosdz

        ds2dx = ds2ddz * ddzetadx + ds2drhos * drhosdx
        ds2dy = ds2ddz * ddzetady + ds2drhos * drhosdy
        ds2dz = ds2ddz * ddzetadz + ds2drhos * drhosdz

        s1ts2 = s1 * s2
        s1ps2 = s1 + s2
        s1ps2sq = s1ps2^2
        fac1 = sqrt(s1ps2sq - (2 * bi)^2)
        as0 = fac1 / (s1ts2 * s1ps2sq)
        term1 = 1 / (s1ts2 * s1ps2 * fac1)
        fac2 = as0 / s1ps2sq
        dasds1 = term1 - fac2 / s1 * (s2^2 + s1 * (3 * s1 + 4 * s2))
        dasds2 = term1 - fac2 / s2 * (s1^2 + s2 * (3 * s2 + 4 * s1))

        dasdx = dasds1 * ds1dx + dasds2 * ds2dx
        dasdy = dasds1 * ds1dy + dasds2 * ds2dy
        dasdz = dasds1 * ds1dz + dasds2 * ds2dz

        bx += F_TD[i] * ((2 * as0 + y * dasdy) * spss - (xs - XSHIFT_TD) * dasdz + as0 * dpsrr * (y^2 * cpss + z * zsww))
        by -= F_TD[i] * y * (as0 * dpsrr * xs + dasdz * cpss + dasdx * spss)
        bz += F_TD[i] * ((2 * as0 + y * dasdy) * cpss + (xs - XSHIFT_TD) * dasdx - as0 * dpsrr * (x * zsww + y^2 * spss))
    end

    return bx, by, bz
end

function tail87(x, z, warp)
    rps = warp.rps
    warpval = warp.warp
    dd = 3.0
    hpi = pi / 2
    rt = 40.0
    xn = -10.0
    x1 = -1.261
    x2 = -0.663
    b0 = 0.391734
    b1 = 5.89715
    b2 = 24.6833
    tscale = 1.0

    b1 *= tscale
    b2 *= tscale^2

    xn21 = (xn - x1)^2
    xnr = 1 / (xn - x2)
    adln = -log(xnr^2 * xn21)

    zs = z - rps + warpval
    zp = z - rt
    zm = z + rt

    xnx = xn - x
    xnx2 = xnx^2
    xc1 = x - x1
    xc2 = x - x2
    xc22 = xc2^2
    xr2 = xc2 * xnr
    xc12 = xc1^2
    d2 = dd^2
    b20 = zs^2 + d2
    b2p = zp^2 + d2
    b2m = zm^2 + d2
    b = sqrt(b20)
    bp = sqrt(b2p)
    bm = sqrt(b2m)
    xa1 = xc12 + b20
    xap1 = xc12 + b2p
    xam1 = xc12 + b2m
    xa2 = 1 / (xc22 + b20)
    xap2 = 1 / (xc22 + b2p)
    xam2 = 1 / (xc22 + b2m)
    xna = xnx2 + b20
    xnap = xnx2 + b2p
    xnam = xnx2 + b2m
    f = b20 - xc22
    fp = b2p - xc22
    fm = b2m - xc22
    xln1 = log(xn21 / xna)
    xlnp1 = log(xn21 / xnap)
    xlnm1 = log(xn21 / xnam)
    xln2 = xln1 + adln
    xlnp2 = xlnp1 + adln
    xlnm2 = xlnm1 + adln
    aln = 0.25 * (xlnp1 + xlnm1 - 2 * xln1)
    s0 = (atan(xnx / b) + hpi) / b
    s0p = (atan(xnx / bp) + hpi) / bp
    s0m = (atan(xnx / bm) + hpi) / bm
    s1 = (xln1 * 0.5 + xc1 * s0) / xa1
    s1p = (xlnp1 * 0.5 + xc1 * s0p) / xap1
    s1m = (xlnm1 * 0.5 + xc1 * s0m) / xam1
    s2 = (xc2 * xa2 * xln2 - xnr - f * xa2 * s0) * xa2
    s2p = (xc2 * xap2 * xlnp2 - xnr - fp * xap2 * s0p) * xap2
    s2m = (xc2 * xam2 * xlnm2 - xnr - fm * xam2 * s0m) * xam2
    g1 = (b20 * s0 - 0.5 * xc1 * xln1) / xa1
    g1p = (b2p * s0p - 0.5 * xc1 * xlnp1) / xap1
    g1m = (b2m * s0m - 0.5 * xc1 * xlnm1) / xam1
    g2 = ((0.5 * f * xln2 + 2 * s0 * b20 * xc2) * xa2 + xr2) * xa2
    g2p = ((0.5 * fp * xlnp2 + 2 * s0p * b2p * xc2) * xap2 + xr2) * xap2
    g2m = ((0.5 * fm * xlnm2 + 2 * s0m * b2m * xc2) * xam2 + xr2) * xam2

    bx = b0 * (zs * s0 - 0.5 * (zp * s0p + zm * s0m)) +
        b1 * (zs * s1 - 0.5 * (zp * s1p + zm * s1m)) +
        b2 * (zs * s2 - 0.5 * (zp * s2p + zm * s2m))
    bz = b0 * aln + b1 * (g1 - 0.5 * (g1p + g1m)) + b2 * (g2 - 0.5 * (g2p + g2m))

    return bx, bz
end


function birk1tot_02(ps, x, y, z)
    rh, dr = 9.0, 4.0
    xltday, xltnght = 78.0, 70.0
    dtet0 = 0.034906
    tnoonn = (90 - xltday) * 0.01745329
    tnoons = pi - tnoonn
    dtetdn = (xltday - xltnght) * 0.01745329

    sps = sin(ps)
    r2 = x^2 + y^2 + z^2
    r = sqrt(r2)
    r3 = r * r2

    rmrh = r - rh
    rprh = r + rh
    dr2 = dr^2
    sqm = sqrt(rmrh^2 + dr2)
    sqp = sqrt(rprh^2 + dr2)
    c = sqp - sqm
    q = sqrt((rh + 1)^2 + dr2) - sqrt((rh - 1)^2 + dr2)
    spsas = sps / r * c / q
    cpsas = sqrt(1 - spsas^2)
    xas = x * cpsas - z * spsas
    zas = x * spsas + z * cpsas
    pas = (xas != 0 || y != 0) ? atan(y, xas) : 0.0
    tas = atan(sqrt(xas^2 + y^2), zas)
    stas = sin(tas)
    f = stas / (stas^6 * (1 - r3) + r3)^(1 / 6)
    tet0 = asin(f)
    if tas > 1.5707963
        tet0 = pi - tet0
    end
    dtet = dtetdn * sin(pas * 0.5)^2
    tetr1n = tnoonn + dtet
    tetr1s = tnoons - dtet

    loc = if tet0 < tetr1n - dtet0 || tet0 > tetr1s + dtet0
        1
    elseif tet0 > tetr1n + dtet0 && tet0 < tetr1s - dtet0
        2
    elseif tet0 >= tetr1n - dtet0 && tet0 <= tetr1n + dtet0
        3
    elseif tet0 >= tetr1s - dtet0 && tet0 <= tetr1s + dtet0
        4
    else
        error("birk1tot_02: invalid region")
    end

    b_diploop1 = if loc in (1, 3, 4)
        @no_escape begin
            d = @alloc(Float64, 3, 26)
            diploop1!(d, x, y, z, ps, rh, dr)
            _bx, _by, _bz = 0.0, 0.0, 0.0
            for i in eachindex(BIRK1_C1)
                _bx += BIRK1_C1[i] * d[1, i]
                _by += BIRK1_C1[i] * d[2, i]
                _bz += BIRK1_C1[i] * d[3, i]
            end
            _bx, _by, _bz
        end
    else
        (0.0, 0.0, 0.0)
    end

    b_condip1 = if loc in (2, 3, 4)
        @no_escape begin
            d2 = @alloc(Float64, 3, 79)
            condip1!(d2, x, y, z, ps)
            _bx, _by, _bz = 0.0, 0.0, 0.0
            for i in eachindex(BIRK1_C2)
                _bx += BIRK1_C2[i] * d2[1, i]
                _by += BIRK1_C2[i] * d2[2, i]
                _bz += BIRK1_C2[i] * d2[3, i]
            end
            _bx, _by, _bz
        end
    else
        (0.0, 0.0, 0.0)
    end

    b = if loc == 1
        b_diploop1
    elseif loc == 2
        b_condip1
    elseif loc == 3
        t01 = tetr1n - dtet0
        t02 = tetr1n + dtet0
        sqr = sqrt(r)
        st01as = sqr / (r3 + 1 / sin(t01)^6 - 1)^(1 / 6)
        st02as = sqr / (r3 + 1 / sin(t02)^6 - 1)^(1 / 6)
        ct01as = sqrt(1 - st01as^2)
        ct02as = sqrt(1 - st02as^2)
        sinpas, cospas = sincos(pas)
        xas1 = r * st01as * cospas
        y1 = r * st01as * sinpas
        zas1 = r * ct01as
        x1 = xas1 * cpsas + zas1 * spsas
        z1 = -xas1 * spsas + zas1 * cpsas
        xas2 = r * st02as * cospas
        y2 = r * st02as * sinpas
        zas2 = r * ct02as
        x2 = xas2 * cpsas + zas2 * spsas
        z2 = -xas2 * spsas + zas2 * cpsas
        ss = sqrt((x2 - x1)^2 + (y2 - y1)^2 + (z2 - z1)^2)
        ds = sqrt((x - x1)^2 + (y - y1)^2 + (z - z1)^2)
        frac = ds / ss
        @. b_diploop1 * (1 - frac) + b_condip1 * frac
    elseif loc == 4
        t01 = tetr1s - dtet0
        t02 = tetr1s + dtet0
        sqr = sqrt(r)
        st01as = sqr / (r3 + 1 / sin(t01)^6 - 1)^(1 / 6)
        st02as = sqr / (r3 + 1 / sin(t02)^6 - 1)^(1 / 6)
        ct01as = -sqrt(1 - st01as^2)
        ct02as = -sqrt(1 - st02as^2)
        sinpas, cospas = sincos(pas)
        xas1 = r * st01as * cospas
        y1 = r * st01as * sinpas
        zas1 = r * ct01as
        x1 = xas1 * cpsas + zas1 * spsas
        z1 = -xas1 * spsas + zas1 * cpsas
        xas2 = r * st02as * cospas
        y2 = r * st02as * sinpas
        zas2 = r * ct02as
        x2 = xas2 * cpsas + zas2 * spsas
        z2 = -xas2 * spsas + zas2 * cpsas
        ss = sqrt((x2 - x1)^2 + (y2 - y1)^2 + (z2 - z1)^2)
        ds = sqrt((x - x1)^2 + (y - y1)^2 + (z - z1)^2)
        frac = ds / ss
        @. b_condip1 * (1 - frac) + b_diploop1 * frac
    end

    return b .+ birk1shld(ps, x, y, z)
end

function diploop1!(d, x, y, z, ps, rh, dr)
    sps = sin(ps)

    for i in 1:12
        r2 = (BIRK1_XX1[i] * BIRK1_DIPX)^2 + (BIRK1_YY1[i] * BIRK1_DIPY)^2
        r = sqrt(r2)
        rmrh = r - rh
        rprh = r + rh
        dr2 = dr^2
        sqm = sqrt(rmrh^2 + dr2)
        sqp = sqrt(rprh^2 + dr2)
        c = sqp - sqm
        q = sqrt((rh + 1)^2 + dr2) - sqrt((rh - 1)^2 + dr2)
        spsas = sps / r * c / q
        cpsas = sqrt(1 - spsas^2)
        xd = (BIRK1_XX1[i] * BIRK1_DIPX) * cpsas
        yd = (BIRK1_YY1[i] * BIRK1_DIPY)
        zd = -(BIRK1_XX1[i] * BIRK1_DIPX) * spsas
        bx1x, by1x, bz1x, bx1y, by1y, bz1y, bx1z, by1z, bz1z = dipxyz(x - xd, y - yd, z - zd)
        if abs(yd) > 1.0e-10
            bx2x, by2x, bz2x, bx2y, by2y, bz2y, bx2z, by2z, bz2z = dipxyz(x - xd, y + yd, z - zd)
        else
            bx2x = by2x = bz2x = 0.0
            bx2y = by2y = bz2y = 0.0
            bx2z = by2z = bz2z = 0.0
        end
        d[1, i] = bx1z + bx2z
        d[2, i] = by1z + by2z
        d[3, i] = bz1z + bz2z
        d[1, i + 12] = (bx1x + bx2x) * sps
        d[2, i + 12] = (by1x + by2x) * sps
        d[3, i + 12] = (bz1x + bz2x) * sps
    end

    r2 = (BIRK1_XCENTRE[1] + BIRK1_RADIUS[1])^2
    r = sqrt(r2)
    rmrh = r - rh
    rprh = r + rh
    dr2 = dr^2
    sqm = sqrt(rmrh^2 + dr2)
    sqp = sqrt(rprh^2 + dr2)
    c = sqp - sqm
    q = sqrt((rh + 1)^2 + dr2) - sqrt((rh - 1)^2 + dr2)
    spsas = sps / r * c / q
    cpsas = sqrt(1 - spsas^2)
    xoct1 = x * cpsas - z * spsas
    yoct1 = y
    zoct1 = x * spsas + z * cpsas
    bxoct1, byoct1, bzoct1 = crosslp(xoct1, yoct1, zoct1, BIRK1_XCENTRE[1], BIRK1_RADIUS[1], BIRK1_TILT)
    d[1, 25] = bxoct1 * cpsas + bzoct1 * spsas
    d[2, 25] = byoct1
    d[3, 25] = -bxoct1 * spsas + bzoct1 * cpsas

    r2 = (BIRK1_RADIUS[2] - BIRK1_XCENTRE[2])^2
    r = sqrt(r2)
    rmrh = r - rh
    rprh = r + rh
    dr2 = dr^2
    sqm = sqrt(rmrh^2 + dr2)
    sqp = sqrt(rprh^2 + dr2)
    c = sqp - sqm
    q = sqrt((rh + 1)^2 + dr2) - sqrt((rh - 1)^2 + dr2)
    spsas = sps / r * c / q
    cpsas = sqrt(1 - spsas^2)
    xoct2 = x * cpsas - z * spsas - BIRK1_XCENTRE[2]
    yoct2 = y
    zoct2 = x * spsas + z * cpsas
    bx, by, bz = circle(xoct2, yoct2, zoct2, BIRK1_RADIUS[2])
    d[1, 26] = bx * cpsas + bz * spsas
    d[2, 26] = by
    d[3, 26] = -bx * spsas + bz * cpsas

    return d
end

function dipxyz(x, y, z)
    x2 = x^2
    y2 = y^2
    z2 = z^2
    r2 = x2 + y2 + z2

    xmr5 = 30574 / (r2 * r2 * sqrt(r2))
    xmr53 = 3 * xmr5
    bxx = xmr5 * (3 * x2 - r2)
    byx = xmr53 * x * y
    bzx = xmr53 * x * z

    bxy = byx
    byy = xmr5 * (3 * y2 - r2)
    bzy = xmr53 * y * z

    bxz = bzx
    byz = bzy
    bzz = xmr5 * (3 * z2 - r2)

    return bxx, byx, bzx, bxy, byy, bzy, bxz, byz, bzz
end

function crosslp(x, y, z, xc, rl, al)
    sal, cal = sincos(al)

    y1 = y * cal - z * sal
    z1 = y * sal + z * cal
    y2 = y * cal + z * sal
    z2 = -y * sal + z * cal
    bx1, by1, bz1 = circle(x - xc, y1, z1, rl)
    bx2, by2, bz2 = circle(x - xc, y2, z2, rl)
    bx = bx1 + bx2
    by = (by1 + by2) * cal + (bz1 - bz2) * sal
    bz = -(by1 - by2) * sal + (bz1 + bz2) * cal

    return bx, by, bz
end

function circle(x, y, z, rl)
    rho2 = x^2 + y^2
    rho = sqrt(rho2)
    r22 = z^2 + (rho + rl)^2
    r2 = sqrt(r22)
    r12 = r22 - 4 * rho * rl
    r32 = 0.5 * (r12 + r22)
    xk2 = 1 - r12 / r22
    xk2s = 1 - xk2
    dl = log(1 / xk2s)
    k = 1.38629436112 + xk2s * (0.09666344259 + xk2s * (0.03590092383 + xk2s * (0.03742563713 + xk2s * 0.01451196212))) +
        dl * (0.5 + xk2s * (0.12498593597 + xk2s * (0.06880248576 + xk2s * (0.03328355346 + xk2s * 0.00441787012))))
    e = 1 + xk2s * (0.44325141463 + xk2s * (0.0626060122 + xk2s * (0.04757383546 + xk2s * 0.01736506451))) +
        dl * xk2s * (0.2499836831 + xk2s * (0.09200180037 + xk2s * (0.04069697526 + xk2s * 0.00526449639)))

    brho = if rho > 1.0e-6
        z / (rho2 * r2) * (r32 / r12 * e - k)
    else
        pi * rl / r2 * (rl - rho) / r12 * z / (r32 - rho2)
    end

    bx = brho * x
    by = brho * y
    bz = (k - e * (r32 - 2 * rl^2) / r12) / r2

    return bx, by, bz
end

function condip1!(d, x, y, z, ps)
    sps, cps = sincos(ps)

    xsm = x * cps - z * sps - BIRK1_DX
    zsm = z * cps + x * sps
    ro2 = xsm^2 + y^2
    ro = sqrt(ro2)

    cf1 = xsm / ro
    sf1 = y / ro
    cf2 = cf1^2 - sf1^2
    sf2 = 2 * sf1 * cf1
    cf3 = cf2 * cf1 - sf2 * sf1
    sf3 = sf2 * cf1 + cf2 * sf1
    cf4 = cf3 * cf1 - sf3 * sf1
    sf4 = sf3 * cf1 + cf3 * sf1
    cf5 = cf4 * cf1 - sf4 * sf1
    sf5 = sf4 * cf1 + cf4 * sf1

    cf = (cf1, cf2, cf3, cf4, cf5)
    sf = (sf1, sf2, sf3, sf4, sf5)

    r2 = ro2 + zsm^2
    r = sqrt(r2)
    c = zsm / r
    s = ro / r
    ch = sqrt(0.5 * (1 + c))
    sh = sqrt(0.5 * (1 - c))
    tnh = sh / ch
    cnh = 1 / tnh

    for m in 1:5
        m1 = m
        bt = m1 * cf[m] / (r * s) * (tnh^m1 + cnh^m1)
        bf = -0.5 * m1 * sf[m] / r * (tnh^(m1 - 1) / ch^2 - cnh^(m1 - 1) / sh^2)
        bxsm = bt * c * cf1 - bf * sf1
        by = bt * c * sf1 + bf * cf1
        bzsm = -bt * s
        d[1, m] = bxsm * cps + bzsm * sps
        d[2, m] = by
        d[3, m] = -bxsm * sps + bzsm * cps
    end

    xsm = x * cps - z * sps
    zsm = z * cps + x * sps

    for i in 1:9
        if i == 3 || i == 5 || i == 6
            xd = BIRK1_XX2[i] * BIRK1_SCALEIN
            yd = BIRK1_YY2[i] * BIRK1_SCALEIN
        else
            xd = BIRK1_XX2[i] * BIRK1_SCALEOUT
            yd = BIRK1_YY2[i] * BIRK1_SCALEOUT
        end
        zd = BIRK1_ZZ2[i]
        bx1x, by1x, bz1x, bx1y, by1y, bz1y, bx1z, by1z, bz1z = dipxyz(xsm - xd, y - yd, zsm - zd)
        bx2x, by2x, bz2x, bx2y, by2y, bz2y, bx2z, by2z, bz2z = dipxyz(xsm - xd, y + yd, zsm - zd)
        bx3x, by3x, bz3x, bx3y, by3y, bz3y, bx3z, by3z, bz3z = dipxyz(xsm - xd, y - yd, zsm + zd)
        bx4x, by4x, bz4x, bx4y, by4y, bz4y, bx4z, by4z, bz4z = dipxyz(xsm - xd, y + yd, zsm + zd)

        ix = i * 3 + 3
        iy = ix + 1
        iz = iy + 1

        d[1, ix] = (bx1x + bx2x - bx3x - bx4x) * cps + (bz1x + bz2x - bz3x - bz4x) * sps
        d[2, ix] = by1x + by2x - by3x - by4x
        d[3, ix] = (bz1x + bz2x - bz3x - bz4x) * cps - (bx1x + bx2x - bx3x - bx4x) * sps

        d[1, iy] = (bx1y - bx2y - bx3y + bx4y) * cps + (bz1y - bz2y - bz3y + bz4y) * sps
        d[2, iy] = by1y - by2y - by3y + by4y
        d[3, iy] = (bz1y - bz2y - bz3y + bz4y) * cps - (bx1y - bx2y - bx3y + bx4y) * sps

        d[1, iz] = (bx1z + bx2z + bx3z + bx4z) * cps + (bz1z + bz2z + bz3z + bz4z) * sps
        d[2, iz] = by1z + by2z + by3z + by4z
        d[3, iz] = (bz1z + bz2z + bz3z + bz4z) * cps - (bx1z + bx2z + bx3z + bx4z) * sps

        ix += 27
        iy += 27
        iz += 27

        d[1, ix] = sps * ((bx1x + bx2x + bx3x + bx4x) * cps + (bz1x + bz2x + bz3x + bz4x) * sps)
        d[2, ix] = sps * (by1x + by2x + by3x + by4x)
        d[3, ix] = sps * ((bz1x + bz2x + bz3x + bz4x) * cps - (bx1x + bx2x + bx3x + bx4x) * sps)

        d[1, iy] = sps * ((bx1y - bx2y + bx3y - bx4y) * cps + (bz1y - bz2y + bz3y - bz4y) * sps)
        d[2, iy] = sps * (by1y - by2y + by3y - by4y)
        d[3, iy] = sps * ((bz1y - bz2y + bz3y - bz4y) * cps - (bx1y - bx2y + bx3y - bx4y) * sps)

        d[1, iz] = sps * ((bx1z + bx2z - bx3z - bx4z) * cps + (bz1z + bz2z - bz3z - bz4z) * sps)
        d[2, iz] = sps * (by1z + by2z - by3z - by4z)
        d[3, iz] = sps * ((bz1z + bz2z - bz3z - bz4z) * cps - (bx1z + bx2z - bx3z - bx4z) * sps)
    end

    for i in 1:5
        zd = BIRK1_ZZ2[i + 9]
        bx1x, by1x, bz1x, bx1y, by1y, bz1y, bx1z, by1z, bz1z = dipxyz(xsm, y, zsm - zd)
        bx2x, by2x, bz2x, bx2y, by2y, bz2y, bx2z, by2z, bz2z = dipxyz(xsm, y, zsm + zd)
        ix = 58 + i * 2
        iz = ix + 1
        d[1, ix] = (bx1x - bx2x) * cps + (bz1x - bz2x) * sps
        d[2, ix] = by1x - by2x
        d[3, ix] = (bz1x - bz2x) * cps - (bx1x - bx2x) * sps

        d[1, iz] = (bx1z + bx2z) * cps + (bz1z + bz2z) * sps
        d[2, iz] = by1z + by2z
        d[3, iz] = (bz1z + bz2z) * cps - (bx1z + bx2z) * sps

        ix += 10
        iz += 10
        d[1, ix] = sps * ((bx1x + bx2x) * cps + (bz1x + bz2x) * sps)
        d[2, ix] = sps * (by1x + by2x)
        d[3, ix] = sps * ((bz1x + bz2x) * cps - (bx1x + bx2x) * sps)

        d[1, iz] = sps * ((bx1z - bx2z) * cps + (bz1z - bz2z) * sps)
        d[2, iz] = sps * (by1z - by2z)
        d[3, iz] = sps * ((bz1z - bz2z) * cps - (bx1z - bx2z) * sps)
    end

    return d
end

@views function birk1shld(ps, x, y, z)
    p1 = NTuple{4}(BIRK1_SHLD[65:68])
    r1 = NTuple{4}(BIRK1_SHLD[69:72])
    q1 = NTuple{4}(BIRK1_SHLD[73:76])
    s1 = NTuple{4}(BIRK1_SHLD[77:80])
    rp = 1.0 ./ p1
    rr = 1.0 ./ r1
    rq = 1.0 ./ q1
    rs = 1.0 ./ s1

    bx = 0.0
    by = 0.0
    bz = 0.0
    sps, cps = sincos(ps)
    s3ps = 4 * cps^2 - 1
    l = 1
    for m in 1:2
        for i in 1:4
            sypi, cypi = sincos(y * rp[i])
            syqi, cyqi = sincos(y * rq[i])
            for k in 1:4
                szrk, czrk = sincos(z * rr[k])
                szsk, czsk = sincos(z * rs[k])
                sqpr = sqrt(rp[i]^2 + rr[k]^2)
                sqqs = sqrt(rq[i]^2 + rs[k]^2)
                epr = exp(x * sqpr)
                eqs = exp(x * sqqs)
                if m == 1
                    hx = -sqpr * epr * cypi * szrk
                    hy = rp[i] * epr * sypi * szrk
                    hz = -rr[k] * epr * cypi * czrk
                    tmp = BIRK1_SHLD[l + 1] * cps + BIRK1_SHLD[l]
                else
                    hx = -sps * sqqs * eqs * cyqi * czsk
                    hy = sps * rq[i] * eqs * syqi * czsk
                    hz = sps * rs[k] * eqs * cyqi * szsk
                    tmp = BIRK1_SHLD[l + 1] * s3ps + BIRK1_SHLD[l]
                end
                bx += hx * tmp
                by += hy * tmp
                bz += hz * tmp
                l += 2
            end
        end
    end

    return bx, by, bz
end

function birk2tot_02(ps, x, y, z)
    return birk2shl(x, y, z, ps) .+ r2_birk(x, y, z, ps)
end

@views function birk2shl(x, y, z, ps)
    p = BIRK2_SHL[17:18]
    r = BIRK2_SHL[19:20]
    q = BIRK2_SHL[21:22]
    s = BIRK2_SHL[23:24]

    sps, cps = sincos(ps)
    s3ps = 4 * cps^2 - 1
    hx = 0.0
    hy = 0.0
    hz = 0.0
    l = 1
    for m in 1:2, i in 1:2
        sypi, cypi = sincos(y / p[i])
        syqi, cyqi = sincos(y / q[i])
        for k in 1:2
            sqpr = sqrt(1 / p[i]^2 + 1 / r[k]^2)
            sqqs = sqrt(1 / q[i]^2 + 1 / s[k]^2)
            epr = exp(x * sqpr)
            eqs = exp(x * sqqs)
            if m == 1
                szrk, czrk = sincos(z / r[k])
                dx = -sqpr * epr * cypi * szrk
                dy = epr / p[i] * sypi * szrk
                dz = -epr / r[k] * cypi * czrk
                tmp = cps
            else
                szsk, czsk = sincos(z / s[k])
                dx = -sps * sqqs * eqs * cyqi * czsk
                dy = sps * eqs / q[i] * syqi * czsk
                dz = sps * eqs / s[k] * cyqi * szsk
                tmp = s3ps
            end
            tmp = BIRK2_SHL[l] + BIRK2_SHL[l + 1] * tmp
            hx += dx * tmp
            hy += dy * tmp
            hz += dz * tmp
            l += 2
        end
    end
    return hx, hy, hz
end

function r2_birk(x, y, z, ps)
    delarg, delarg1 = 0.03, 0.015
    sps, cps = sincos(ps)
    xsm = x * cps - z * sps
    zsm = z * cps + x * sps
    xks = xksi(xsm, y, zsm)

    bxsm, by, bzsm = if xks < -(delarg + delarg1)
        r2outer(xsm, y, zsm) .* -0.02
    elseif xks < -delarg + delarg1
        f2 = -0.02 * tksi(xks, -delarg, delarg1)
        f1 = -0.02 - f2
        r2outer(xsm, y, zsm) .* f1 .+ r2sheet(xsm, y, zsm) .* f2
    elseif xks < delarg - delarg1
        r2sheet(xsm, y, zsm) .* -0.02
    elseif xks < delarg + delarg1
        f1 = -0.02 * tksi(xks, delarg, delarg1)
        f2 = -0.02 - f1
        r2inner(xsm, y, zsm) .* f1 .+ r2sheet(xsm, y, zsm) .* f2
    else
        r2inner(xsm, y, zsm) .* -0.02
    end

    bx = bxsm * cps + bzsm * sps
    bz = bzsm * cps - bxsm * sps
    return bx, by, bz
end

function xksi(x, y, z)
    a11a12, a21a22, a41a42, a51a52, a61a62, b11b12, b21b22, c61c62, c71c72, r0, dr =
        0.305662, -0.383593, 0.2677733, -0.097656, -0.636034, -0.359862,
        0.424706, -0.126366, 0.292578, 1.21563, 7.50937
    tnoon, dteta = 0.3665191, 0.09599309

    dr2 = dr^2
    r2 = x^2 + y^2 + z^2
    r = sqrt(r2)
    xr = x / r
    yr = y / r
    zr = z / r

    pr = r < r0 ? 0.0 : sqrt((r - r0)^2 + dr2) - dr
    f = x + pr * (a11a12 + a21a22 * xr + a41a42 * xr^2 + a51a52 * yr^2 + a61a62 * zr^2)
    g = y + pr * (b11b12 * yr + b21b22 * xr * yr)
    h = z + pr * (c61c62 * zr + c71c72 * xr * zr)
    g2 = g^2

    fgh = f^2 + g2 + h^2
    fgh32 = sqrt(fgh)^3
    fchsg2 = f^2 + g2

    if fchsg2 < 1.0e-5
        return -1.0
    else
        sqfchsg2 = sqrt(fchsg2)
        alpha = fchsg2 / fgh32
        theta = tnoon + 0.5 * dteta * (1 - f / sqfchsg2)
        phi = sin(theta)^2
        return alpha - phi
    end
end

function tksi(xksi, xks0, dxksi)
    tdz3 = 2 * dxksi^3
    return if xksi - xks0 < -dxksi
        0.0
    elseif xksi < xks0
        br3 = (xksi - xks0 + dxksi)^3
        1.5 * br3 / (tdz3 + br3)
    elseif xksi - xks0 < dxksi
        br3 = (xksi - xks0 - dxksi)^3
        1.0 + 1.5 * br3 / (tdz3 - br3)
    else
        1.0
    end
end

fexp(s, a) = a < 0 ? sqrt(-2 * a * MathConstants.e) * s * exp(a * s^2) : s * exp(a * (s^2 - 1))
fexp1(s, a) = a <= 0 ? exp(a * s^2) : exp(a * (s^2 - 1))

function r2outer(x, y, z)
    pn = R2_OUTER_PN
    dbx1, dby1, dbz1 = crosslp(x, y, z, pn[1], pn[2], pn[3])
    dbx2, dby2, dbz2 = crosslp(x, y, z, pn[4], pn[5], pn[6])
    dbx3, dby3, dbz3 = crosslp(x, y, z, pn[7], pn[8], pn[9])
    dbx4, dby4, dbz4 = circle(x - pn[10], y, z, pn[11])
    dbx5, dby5, dbz5 = loops4(x, y, z, pn[12], pn[13], pn[14], pn[15], pn[16], pn[17])

    bx = R2_OUTER_PL[1] * dbx1 + R2_OUTER_PL[2] * dbx2 + R2_OUTER_PL[3] * dbx3 + R2_OUTER_PL[4] * dbx4 + R2_OUTER_PL[5] * dbx5
    by = R2_OUTER_PL[1] * dby1 + R2_OUTER_PL[2] * dby2 + R2_OUTER_PL[3] * dby3 + R2_OUTER_PL[4] * dby4 + R2_OUTER_PL[5] * dby5
    bz = R2_OUTER_PL[1] * dbz1 + R2_OUTER_PL[2] * dbz2 + R2_OUTER_PL[3] * dbz3 + R2_OUTER_PL[4] * dbz4 + R2_OUTER_PL[5] * dbz5
    return bx, by, bz
end

function loops4(x, y, z, xc, yc, zc, r, theta, phi)
    st, ct = sincos(theta)
    sp, cp = sincos(phi)

    xs = (x - xc) * cp + (y - yc) * sp
    yss = (y - yc) * cp - (x - xc) * sp
    zs = z - zc
    xss = xs * ct - zs * st
    zss = zs * ct + xs * st
    bxss, bys, bzss = circle(xss, yss, zss, r)
    bxs = bxss * ct + bzss * st
    bz1 = bzss * ct - bxss * st
    bx1 = bxs * cp - bys * sp
    by1 = bxs * sp + bys * cp

    xs = (x - xc) * cp - (y + yc) * sp
    yss = (y + yc) * cp + (x - xc) * sp
    zs = z - zc
    xss = xs * ct - zs * st
    zss = zs * ct + xs * st
    bxss, bys, bzss = circle(xss, yss, zss, r)
    bxs = bxss * ct + bzss * st
    bz2 = bzss * ct - bxss * st
    bx2 = bxs * cp + bys * sp
    by2 = -bxs * sp + bys * cp

    xs = -(x - xc) * cp + (y + yc) * sp
    yss = -(y + yc) * cp - (x - xc) * sp
    zs = z + zc
    xss = xs * ct - zs * st
    zss = zs * ct + xs * st
    bxss, bys, bzss = circle(xss, yss, zss, r)
    bxs = bxss * ct + bzss * st
    bz3 = bzss * ct - bxss * st
    bx3 = -bxs * cp - bys * sp
    by3 = bxs * sp - bys * cp

    xs = -(x - xc) * cp - (y - yc) * sp
    yss = -(y - yc) * cp + (x - xc) * sp
    zs = z + zc
    xss = xs * ct - zs * st
    zss = zs * ct + xs * st
    bxss, bys, bzss = circle(xss, yss, zss, r)
    bxs = bxss * ct + bzss * st
    bz4 = bzss * ct - bxss * st
    bx4 = -bxs * cp + bys * sp
    by4 = -bxs * sp - bys * cp

    bx = bx1 + bx2 + bx3 + bx4
    by = by1 + by2 + by3 + by4
    bz = bz1 + bz2 + bz3 + bz4

    return bx, by, bz
end

function r2sheet(x, y, z)
    xks = xksi(x, y, z)
    t1x = xks / sqrt(xks^2 + R2_SHEET_PNONX[6]^2)
    t2x = R2_SHEET_PNONX[7]^3 / sqrt(xks^2 + R2_SHEET_PNONX[7]^2)^3
    t3x = xks / sqrt(xks^2 + R2_SHEET_PNONX[8]^2)^5 * 3.493856 * R2_SHEET_PNONX[8]^4

    t1y = xks / sqrt(xks^2 + R2_SHEET_PNONY[6]^2)
    t2y = R2_SHEET_PNONY[7]^3 / sqrt(xks^2 + R2_SHEET_PNONY[7]^2)^3
    t3y = xks / sqrt(xks^2 + R2_SHEET_PNONY[8]^2)^5 * 3.493856 * R2_SHEET_PNONY[8]^4

    t1z = xks / sqrt(xks^2 + R2_SHEET_PNONZ[6]^2)
    t2z = R2_SHEET_PNONZ[7]^3 / sqrt(xks^2 + R2_SHEET_PNONZ[7]^2)^3
    t3z = xks / sqrt(xks^2 + R2_SHEET_PNONZ[8]^2)^5 * 3.493856 * R2_SHEET_PNONZ[8]^4

    rho2 = x^2 + y^2
    r = sqrt(rho2 + z^2)
    rho = sqrt(rho2)
    if rho < 1.0e-8
        c1p = 1.0
        s1p = 0.0
    else
        c1p = x / rho
        s1p = y / rho
    end
    s2p = 2 * s1p * c1p
    c2p = c1p^2 - s1p^2
    s3p = s2p * c1p + c2p * s1p
    c3p = c2p * c1p - s2p * s1p
    s4p = s3p * c1p + c3p * s1p
    ct = z / r

    s1 = fexp(ct, R2_SHEET_PNONX[1])
    s2 = fexp(ct, R2_SHEET_PNONX[2])
    s3 = fexp(ct, R2_SHEET_PNONX[3])
    s4 = fexp(ct, R2_SHEET_PNONX[4])
    s5 = fexp(ct, R2_SHEET_PNONX[5])
    bx = s1 * (
        (R2_SHEET_A[1] + R2_SHEET_A[2] * t1x + R2_SHEET_A[3] * t2x + R2_SHEET_A[4] * t3x) +
            c1p * (R2_SHEET_A[5] + R2_SHEET_A[6] * t1x + R2_SHEET_A[7] * t2x + R2_SHEET_A[8] * t3x) +
            c2p * (R2_SHEET_A[9] + R2_SHEET_A[10] * t1x + R2_SHEET_A[11] * t2x + R2_SHEET_A[12] * t3x) +
            c3p * (R2_SHEET_A[13] + R2_SHEET_A[14] * t1x + R2_SHEET_A[15] * t2x + R2_SHEET_A[16] * t3x)
    ) +
        s2 * (
        (R2_SHEET_A[17] + R2_SHEET_A[18] * t1x + R2_SHEET_A[19] * t2x + R2_SHEET_A[20] * t3x) +
            c1p * (R2_SHEET_A[21] + R2_SHEET_A[22] * t1x + R2_SHEET_A[23] * t2x + R2_SHEET_A[24] * t3x) +
            c2p * (R2_SHEET_A[25] + R2_SHEET_A[26] * t1x + R2_SHEET_A[27] * t2x + R2_SHEET_A[28] * t3x) +
            c3p * (R2_SHEET_A[29] + R2_SHEET_A[30] * t1x + R2_SHEET_A[31] * t2x + R2_SHEET_A[32] * t3x)
    ) +
        s3 * (
        (R2_SHEET_A[33] + R2_SHEET_A[34] * t1x + R2_SHEET_A[35] * t2x + R2_SHEET_A[36] * t3x) +
            c1p * (R2_SHEET_A[37] + R2_SHEET_A[38] * t1x + R2_SHEET_A[39] * t2x + R2_SHEET_A[40] * t3x) +
            c2p * (R2_SHEET_A[41] + R2_SHEET_A[42] * t1x + R2_SHEET_A[43] * t2x + R2_SHEET_A[44] * t3x) +
            c3p * (R2_SHEET_A[45] + R2_SHEET_A[46] * t1x + R2_SHEET_A[47] * t2x + R2_SHEET_A[48] * t3x)
    ) +
        s4 * (
        (R2_SHEET_A[49] + R2_SHEET_A[50] * t1x + R2_SHEET_A[51] * t2x + R2_SHEET_A[52] * t3x) +
            c1p * (R2_SHEET_A[53] + R2_SHEET_A[54] * t1x + R2_SHEET_A[55] * t2x + R2_SHEET_A[56] * t3x) +
            c2p * (R2_SHEET_A[57] + R2_SHEET_A[58] * t1x + R2_SHEET_A[59] * t2x + R2_SHEET_A[60] * t3x) +
            c3p * (R2_SHEET_A[61] + R2_SHEET_A[62] * t1x + R2_SHEET_A[63] * t2x + R2_SHEET_A[64] * t3x)
    ) +
        s5 * (
        (R2_SHEET_A[65] + R2_SHEET_A[66] * t1x + R2_SHEET_A[67] * t2x + R2_SHEET_A[68] * t3x) +
            c1p * (R2_SHEET_A[69] + R2_SHEET_A[70] * t1x + R2_SHEET_A[71] * t2x + R2_SHEET_A[72] * t3x) +
            c2p * (R2_SHEET_A[73] + R2_SHEET_A[74] * t1x + R2_SHEET_A[75] * t2x + R2_SHEET_A[76] * t3x) +
            c3p * (R2_SHEET_A[77] + R2_SHEET_A[78] * t1x + R2_SHEET_A[79] * t2x + R2_SHEET_A[80] * t3x)
    )

    s1 = fexp(ct, R2_SHEET_PNONY[1])
    s2 = fexp(ct, R2_SHEET_PNONY[2])
    s3 = fexp(ct, R2_SHEET_PNONY[3])
    s4 = fexp(ct, R2_SHEET_PNONY[4])
    s5 = fexp(ct, R2_SHEET_PNONY[5])
    by = s1 * (
        s1p * (R2_SHEET_B[1] + R2_SHEET_B[2] * t1y + R2_SHEET_B[3] * t2y + R2_SHEET_B[4] * t3y) +
            s2p * (R2_SHEET_B[5] + R2_SHEET_B[6] * t1y + R2_SHEET_B[7] * t2y + R2_SHEET_B[8] * t3y) +
            s3p * (R2_SHEET_B[9] + R2_SHEET_B[10] * t1y + R2_SHEET_B[11] * t2y + R2_SHEET_B[12] * t3y) +
            s4p * (R2_SHEET_B[13] + R2_SHEET_B[14] * t1y + R2_SHEET_B[15] * t2y + R2_SHEET_B[16] * t3y)
    ) +
        s2 * (
        s1p * (R2_SHEET_B[17] + R2_SHEET_B[18] * t1y + R2_SHEET_B[19] * t2y + R2_SHEET_B[20] * t3y) +
            s2p * (R2_SHEET_B[21] + R2_SHEET_B[22] * t1y + R2_SHEET_B[23] * t2y + R2_SHEET_B[24] * t3y) +
            s3p * (R2_SHEET_B[25] + R2_SHEET_B[26] * t1y + R2_SHEET_B[27] * t2y + R2_SHEET_B[28] * t3y) +
            s4p * (R2_SHEET_B[29] + R2_SHEET_B[30] * t1y + R2_SHEET_B[31] * t2y + R2_SHEET_B[32] * t3y)
    ) +
        s3 * (
        s1p * (R2_SHEET_B[33] + R2_SHEET_B[34] * t1y + R2_SHEET_B[35] * t2y + R2_SHEET_B[36] * t3y) +
            s2p * (R2_SHEET_B[37] + R2_SHEET_B[38] * t1y + R2_SHEET_B[39] * t2y + R2_SHEET_B[40] * t3y) +
            s3p * (R2_SHEET_B[41] + R2_SHEET_B[42] * t1y + R2_SHEET_B[43] * t2y + R2_SHEET_B[44] * t3y) +
            s4p * (R2_SHEET_B[45] + R2_SHEET_B[46] * t1y + R2_SHEET_B[47] * t2y + R2_SHEET_B[48] * t3y)
    ) +
        s4 * (
        s1p * (R2_SHEET_B[49] + R2_SHEET_B[50] * t1y + R2_SHEET_B[51] * t2y + R2_SHEET_B[52] * t3y) +
            s2p * (R2_SHEET_B[53] + R2_SHEET_B[54] * t1y + R2_SHEET_B[55] * t2y + R2_SHEET_B[56] * t3y) +
            s3p * (R2_SHEET_B[57] + R2_SHEET_B[58] * t1y + R2_SHEET_B[59] * t2y + R2_SHEET_B[60] * t3y) +
            s4p * (R2_SHEET_B[61] + R2_SHEET_B[62] * t1y + R2_SHEET_B[63] * t2y + R2_SHEET_B[64] * t3y)
    ) +
        s5 * (
        s1p * (R2_SHEET_B[65] + R2_SHEET_B[66] * t1y + R2_SHEET_B[67] * t2y + R2_SHEET_B[68] * t3y) +
            s2p * (R2_SHEET_B[69] + R2_SHEET_B[70] * t1y + R2_SHEET_B[71] * t2y + R2_SHEET_B[72] * t3y) +
            s3p * (R2_SHEET_B[73] + R2_SHEET_B[74] * t1y + R2_SHEET_B[75] * t2y + R2_SHEET_B[76] * t3y) +
            s4p * (R2_SHEET_B[77] + R2_SHEET_B[78] * t1y + R2_SHEET_B[79] * t2y + R2_SHEET_B[80] * t3y)
    )

    s1 = fexp1(ct, R2_SHEET_PNONZ[1])
    s2 = fexp1(ct, R2_SHEET_PNONZ[2])
    s3 = fexp1(ct, R2_SHEET_PNONZ[3])
    s4 = fexp1(ct, R2_SHEET_PNONZ[4])
    s5 = fexp1(ct, R2_SHEET_PNONZ[5])
    bz = s1 * (
        (R2_SHEET_C[1] + R2_SHEET_C[2] * t1z + R2_SHEET_C[3] * t2z + R2_SHEET_C[4] * t3z) +
            c1p * (R2_SHEET_C[5] + R2_SHEET_C[6] * t1z + R2_SHEET_C[7] * t2z + R2_SHEET_C[8] * t3z) +
            c2p * (R2_SHEET_C[9] + R2_SHEET_C[10] * t1z + R2_SHEET_C[11] * t2z + R2_SHEET_C[12] * t3z) +
            c3p * (R2_SHEET_C[13] + R2_SHEET_C[14] * t1z + R2_SHEET_C[15] * t2z + R2_SHEET_C[16] * t3z)
    ) +
        s2 * (
        (R2_SHEET_C[17] + R2_SHEET_C[18] * t1z + R2_SHEET_C[19] * t2z + R2_SHEET_C[20] * t3z) +
            c1p * (R2_SHEET_C[21] + R2_SHEET_C[22] * t1z + R2_SHEET_C[23] * t2z + R2_SHEET_C[24] * t3z) +
            c2p * (R2_SHEET_C[25] + R2_SHEET_C[26] * t1z + R2_SHEET_C[27] * t2z + R2_SHEET_C[28] * t3z) +
            c3p * (R2_SHEET_C[29] + R2_SHEET_C[30] * t1z + R2_SHEET_C[31] * t2z + R2_SHEET_C[32] * t3z)
    ) +
        s3 * (
        (R2_SHEET_C[33] + R2_SHEET_C[34] * t1z + R2_SHEET_C[35] * t2z + R2_SHEET_C[36] * t3z) +
            c1p * (R2_SHEET_C[37] + R2_SHEET_C[38] * t1z + R2_SHEET_C[39] * t2z + R2_SHEET_C[40] * t3z) +
            c2p * (R2_SHEET_C[41] + R2_SHEET_C[42] * t1z + R2_SHEET_C[43] * t2z + R2_SHEET_C[44] * t3z) +
            c3p * (R2_SHEET_C[45] + R2_SHEET_C[46] * t1z + R2_SHEET_C[47] * t2z + R2_SHEET_C[48] * t3z)
    ) +
        s4 * (
        (R2_SHEET_C[49] + R2_SHEET_C[50] * t1z + R2_SHEET_C[51] * t2z + R2_SHEET_C[52] * t3z) +
            c1p * (R2_SHEET_C[53] + R2_SHEET_C[54] * t1z + R2_SHEET_C[55] * t2z + R2_SHEET_C[56] * t3z) +
            c2p * (R2_SHEET_C[57] + R2_SHEET_C[58] * t1z + R2_SHEET_C[59] * t2z + R2_SHEET_C[60] * t3z) +
            c3p * (R2_SHEET_C[61] + R2_SHEET_C[62] * t1z + R2_SHEET_C[63] * t2z + R2_SHEET_C[64] * t3z)
    ) +
        s5 * (
        (R2_SHEET_C[65] + R2_SHEET_C[66] * t1z + R2_SHEET_C[67] * t2z + R2_SHEET_C[68] * t3z) +
            c1p * (R2_SHEET_C[69] + R2_SHEET_C[70] * t1z + R2_SHEET_C[71] * t2z + R2_SHEET_C[72] * t3z) +
            c2p * (R2_SHEET_C[73] + R2_SHEET_C[74] * t1z + R2_SHEET_C[75] * t2z + R2_SHEET_C[76] * t3z) +
            c3p * (R2_SHEET_C[77] + R2_SHEET_C[78] * t1z + R2_SHEET_C[79] * t2z + R2_SHEET_C[80] * t3z)
    )

    return bx, by, bz
end

function r2inner(x, y, z)
    pl = R2_INNER_PL
    pn = R2_INNER_PN
    nmax = 5
    b = @no_escape begin
        cbx = @alloc(Float64, nmax)
        cby = @alloc(Float64, nmax)
        cbz = @alloc(Float64, nmax)
        bconic!(cbx, cby, cbz, x, y, z)
        bx = pl[1] * cbx[1] + pl[2] * cbx[2] + pl[3] * cbx[3] + pl[4] * cbx[4] + pl[5] * cbx[5]
        by = pl[1] * cby[1] + pl[2] * cby[2] + pl[3] * cby[3] + pl[4] * cby[4] + pl[5] * cby[5]
        bz = pl[1] * cbz[1] + pl[2] * cbz[2] + pl[3] * cbz[3] + pl[4] * cbz[4] + pl[5] * cbz[5]
        (bx, by, bz)
    end
    db8 = loops4(x, y, z, pn[1], pn[2], pn[3], pn[4], pn[5], pn[6])
    db6 = dipdistr(x - pn[7], y, z, 0)
    db7 = dipdistr(x - pn[8], y, z, 1)

    return b .+ pl[6] .* db6 .+ pl[7] .* db7 .+ pl[8] .* db8
end

function bconic!(cbx, cby, cbz, x, y, z)
    ro2 = x^2 + y^2
    ro = sqrt(ro2)
    cf = x / ro
    sf = y / ro
    cfm1 = 1.0
    sfm1 = 0.0
    r2 = ro2 + z^2
    r = sqrt(r2)
    c = z / r
    s = ro / r
    ch = sqrt(0.5 * (1 + c))
    sh = sqrt(0.5 * (1 - c))
    tnhm1 = 1.0
    cnhm1 = 1.0
    tnh = sh / ch
    cnh = 1 / tnh

    for m in eachindex(cbx)
        m1 = m
        cfm = cfm1 * cf - sfm1 * sf
        sfm = cfm1 * sf + sfm1 * cf
        cfm1 = cfm
        sfm1 = sfm
        tnhm = tnhm1 * tnh
        cnhm = cnhm1 * cnh
        bt = m1 * cfm / (r * s) * (tnhm + cnhm)
        bf = -0.5 * m1 * sfm / r * (tnhm1 / ch^2 - cnhm1 / sh^2)
        tnhm1 = tnhm
        cnhm1 = cnhm
        cbx[m] = bt * c * cf - bf * sf
        cby[m] = bt * c * sf + bf * cf
        cbz[m] = -bt * s
    end

    return cbx, cby, cbz
end
