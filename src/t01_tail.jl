# T01 Tail Field Functions

# Coefficients for unwarped mode 1
const TAIL_A1 = (
    -25.45869857, 57.3589908, 317.5501869, -2.626756717, -93.38053698,
    -199.6467926, -858.8129729, 34.09192395, 845.4214929, -29.07463068,
    47.10678547, -128.9797943, -781.7512093, 6.165038619, 167.8905046,
    492.068041, 1654.724031, -46.7733792, -1635.922669, 40.86186772,
    -0.1349775602, -0.9661991179e-1, -0.1662302354, 0.002810467517, 0.2487355077,
    0.1025565237, -14.41750229, -0.8185333989, 11.07693629, 0.7569503173,
    -9.655264745, 112.2446542, 777.5948964, -5.745008536, -83.03921993,
    -490.2278695, -1155.004209, 39.0802332, 1172.780574, -39.44349797,
    -14.07211198, -40.41201127, -313.2277343, 2.203920979, 8.232835341,
    197.7065115, 391.2733948, -18.57424451, -437.2779053, 23.04976898,
    11.75673963, 13.60497313, 4.69192706, 18.20923547, 27.59044809,
    6.677425469, 1.398283308, 2.839005878, 31.24817706, 24.53577264,
)

# Coefficients for unwarped mode 2
const TAIL_A2 = (
    -287187.1962, 4970.499233, 410490.1952, -1347.839052, -386370.324,
    3317.98375, -143462.3895, 5706.513767, 171176.2904, 250.888275,
    -506570.8891, 5733.592632, 397975.5842, 9771.762168, -941834.2436,
    7990.97526, 54313.10318, 447.538806, 528046.3449, 12751.04453,
    -21920.98301, -21.05075617, 31971.07875, 3012.641612, -301822.9103,
    -3601.107387, 1797.577552, -6.315855803, 142578.8406, 13161.9364,
    804184.841, -14168.99698, -851926.636, -1890.885671, 972475.6869,
    -8571.862853, 26432.49197, -2554.752298, -482308.3431, -4391.473324,
    105155.916, -1134.62205, -74353.53091, -5382.670711, 695055.0788,
    -916.3365144, -12111.06667, 67.20923358, -367200.9285, -21414.14421,
    14.75567902, 20.7563819, 59.78601609, 16.86431444, 32.58482365,
    23.69472951, 17.24977936, 13.64902647, 68.40989058, 11.67828167,
)

# Taildisk coefficients (from Python)
const TAIL_F = (-71.09346626, -1014.308601, -1272.939359, -3224.935936, -44546.86232)
const TAIL_B = (10.90101242, 12.68393898, 13.51791954, 14.86775017, 15.12306404)
const TAIL_C = (0.7954069972, 0.6716601849, 1.174866319, 2.56524992, 10.0198679)

function deformed(ps, x, y, z)
    rh2, ieps = -5.2, 3
    sps = sin(ps); r2 = x^2 + y^2 + z^2; r = sqrt(r2); zr = z / r
    rh = STATE[].rh0 + rh2 * zr^2
    drhdr = -zr / r * 2.0 * rh2 * zr; drhdz = 2.0 * rh2 * zr / r
    rrh = r / rh; f = 1.0 / (1.0 + rrh^ieps)^(1.0 / ieps)
    dfdr = -rrh^(ieps - 1) * f^(ieps + 1) / rh; dfdrh = -rrh * dfdr
    spsas = sps * f; cpsas = sqrt(1.0 - spsas^2)
    xas, zas = x * cpsas - z * spsas, x * spsas + z * cpsas
    facps = sps / cpsas * (dfdr + dfdrh * drhdr) / r
    psasx, psasy, psasz = facps * x, facps * y, facps * z + sps / cpsas * dfdrh * drhdz
    dxasdx, dxasdy, dxasdz = cpsas - zas * psasx, -zas * psasy, -spsas - zas * psasz
    dzasdx, dzasdy, dzasdz = spsas + xas * psasx, xas * psasy, cpsas + xas * psasz
    fac1 = dxasdz * dzasdy - dxasdy * dzasdz
    fac2 = dxasdx * dzasdz - dxasdz * dzasdx
    fac3 = dzasdx * dxasdy - dxasdx * dzasdy
    bxas1, byas1, bzas1, bxas2, byas2, bzas2 = warped(ps, xas, y, zas)
    bx1 = bxas1 * dzasdz - bzas1 * dxasdz + byas1 * fac1
    by1 = byas1 * fac2
    bz1 = bzas1 * dxasdx - bxas1 * dzasdx + byas1 * fac3
    bx2 = bxas2 * dzasdz - bzas2 * dxasdz + byas2 * fac1
    by2 = byas2 * fac2
    bz2 = bzas2 * dxasdx - bxas2 * dzasdx + byas2 * fac3
    return bx1, by1, bz1, bx2, by2, bz2
end

function warped(ps, x, y, z)
    xl = 20.0; sps = sin(ps); rho2 = y^2 + z^2; rho = sqrt(rho2)
    if y == 0.0 && z == 0.0
        phi, cphi, sphi = 0.0, 1.0, 0.0
    else
        phi = atan(z, y); cphi, sphi = y / rho, z / rho
    end
    rr4l4 = rho / (rho2^2 + xl^4)
    f = phi + STATE[].g * rho2 * rr4l4 * cphi * sps
    dfdphi = 1.0 - STATE[].g * rho2 * rr4l4 * sphi * sps
    dfdrho = STATE[].g * rr4l4^2 * (3.0 * xl^4 - rho2^2) * cphi * sps
    dfdx = 0.0
    cf, sf = cos(f), sin(f); yas, zas = rho * cf, rho * sf
    bx_as1, by_as1, bz_as1, bx_as2, by_as2, bz_as2 = unwarped(x, yas, zas)
    brho_as = by_as1 * cf + bz_as1 * sf; bphi_as = -by_as1 * sf + bz_as1 * cf
    brho_s = brho_as * dfdphi; bphi_s = bphi_as - rho * (bx_as1 * dfdx + brho_as * dfdrho)
    bx1, by1, bz1 = bx_as1 * dfdphi, brho_s * cphi - bphi_s * sphi, brho_s * sphi + bphi_s * cphi
    brho_as = by_as2 * cf + bz_as2 * sf; bphi_as = -by_as2 * sf + bz_as2 * cf
    brho_s = brho_as * dfdphi; bphi_s = bphi_as - rho * (bx_as2 * dfdx + brho_as * dfdrho)
    bx2, by2, bz2 = bx_as2 * dfdphi, brho_s * cphi - bphi_s * sphi, brho_s * sphi + bphi_s * cphi
    return bx1, by1, bz1, bx2, by2, bz2
end

function unwarped(x, y, z)
    # Mode 1 parameters
    deltadx1, alpha1, xshift1 = 1.0, 1.1, 6.0
    # Mode 2 parameters
    deltadx2, alpha2, xshift2 = 0.0, 0.25, 4.0
    xm1, xm2 = -12.0, -12.0

    # Mode 1
    xsc1 = (x - xshift1 - STATE[].dxshift1) * alpha1 - xm1 * (alpha1 - 1.0)
    ysc1 = y * alpha1
    zsc1 = z * alpha1
    d0sc1 = STATE[].d * alpha1

    fx1, fy1, fz1 = taildisk(d0sc1, deltadx1, STATE[].deltady, xsc1, ysc1, zsc1)
    hx1, hy1, hz1 = shlcar5x5(TAIL_A1, x, y, z, STATE[].dxshift1)
    bx1 = fx1 + hx1
    by1 = fy1 + hy1
    bz1 = fz1 + hz1

    # Mode 2
    xsc2 = (x - xshift2 - STATE[].dxshift2) * alpha2 - xm2 * (alpha2 - 1.0)
    ysc2 = y * alpha2
    zsc2 = z * alpha2
    d0sc2 = STATE[].d * alpha2

    fx2, fy2, fz2 = taildisk(d0sc2, deltadx2, STATE[].deltady, xsc2, ysc2, zsc2)
    hx2, hy2, hz2 = shlcar5x5(TAIL_A2, x, y, z, STATE[].dxshift2)
    bx2 = fx2 + hx2
    by2 = fy2 + hy2
    bz2 = fz2 + hz2

    return bx1, by1, bz1, bx2, by2, bz2
end

function taildisk(d0, deltadx, deltady, x, y, z)
    rho = sqrt(x^2 + y^2)
    rho < 1.0e-8 && return 0.0, 0.0, 0.0
    drhodx, drhody = x / rho, y / rho

    dex = exp(x / 7.0)
    d = d0 + deltady * (y / 20.0)^2 + deltadx * dex
    dddy = deltady * y * 0.005
    dddx = deltadx / 7.0 * dex

    dzeta = sqrt(z^2 + d^2)
    ddzetadx = d * dddx / dzeta
    ddzetady = d * dddy / dzeta
    ddzetadz = z / dzeta

    dbx, dby, dbz = 0.0, 0.0, 0.0
    for i in 1:5
        bi, ci = TAIL_B[i], TAIL_C[i]
        s1 = sqrt((rho + bi)^2 + (dzeta + ci)^2)
        s2 = sqrt((rho - bi)^2 + (dzeta + ci)^2)
        ds1drho, ds2drho = (rho + bi) / s1, (rho - bi) / s2
        ds1ddz, ds2ddz = (dzeta + ci) / s1, (dzeta + ci) / s2
        ds1dx = ds1drho * drhodx + ds1ddz * ddzetadx
        ds1dy = ds1drho * drhody + ds1ddz * ddzetady
        ds1dz = ds1ddz * ddzetadz
        ds2dx = ds2drho * drhodx + ds2ddz * ddzetadx
        ds2dy = ds2drho * drhody + ds2ddz * ddzetady
        ds2dz = ds2ddz * ddzetadz
        s1ts2, s1ps2 = s1 * s2, s1 + s2
        s1ps2sq = s1ps2^2
        fac1 = sqrt(s1ps2sq - (2.0 * bi)^2)
        asas = fac1 / (s1ts2 * s1ps2sq)
        dasds1 = (1.0 / (fac1 * s2) - asas / s1ps2 * (s2^2 + s1 * (3.0 * s1 + 4.0 * s2))) / (s1 * s1ps2)
        dasds2 = (1.0 / (fac1 * s1) - asas / s1ps2 * (s1^2 + s2 * (3.0 * s2 + 4.0 * s1))) / (s2 * s1ps2)
        dasdx = dasds1 * ds1dx + dasds2 * ds2dx
        dasdy = dasds1 * ds1dy + dasds2 * ds2dy
        dasdz = dasds1 * ds1dz + dasds2 * ds2dz
        dbx -= TAIL_F[i] * x * dasdz
        dby -= TAIL_F[i] * y * dasdz
        dbz += TAIL_F[i] * (2.0 * asas + x * dasdx + y * dasdy)
    end
    return dbx, dby, dbz
end

function shlcar5x5(a, x, y, z, dshift)
    dhx, dhy, dhz = 0.0, 0.0, 0.0
    l = 0
    for i in 1:5
        rp = 1.0 / a[50 + i]
        cypi, sypi = cos(y * rp), sin(y * rp)
        for k in 1:5
            rr = 1.0 / a[55 + k]
            szrk, czrk = sin(z * rr), cos(z * rr)
            sqpr = sqrt(rp^2 + rr^2)
            epr = exp(x * sqpr)
            dbx = -sqpr * epr * cypi * szrk
            dby = rp * epr * sypi * szrk
            dbz = -rr * epr * cypi * czrk
            l += 1
            coef = a[2 * l - 1] + a[2 * l] * dshift
            dhx += coef * dbx
            dhy += coef * dby
            dhz += coef * dbz
        end
    end
    return dhx, dhy, dhz
end
