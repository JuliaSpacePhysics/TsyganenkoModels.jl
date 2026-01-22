# T01 Core Functions
using ..TsyganenkoModels: dipole

mutable struct T01State
    dxshift1::Float64; dxshift2::Float64; d::Float64; deltady::Float64
    xkappa1::Float64; xkappa2::Float64; sc_sy::Float64; sc_pr::Float64
    phi::Float64; g::Float64; rh0::Float64
end
const STATE = Ref(T01State(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.0))

# Python a[i] -> Julia a[i+1] (0-indexed to 1-indexed)
function t01_extall(a, pdyn, dst, byimf, bzimf, g1, g2, ps, x, y, z)
    rh2 = -5.2
    # a[38] in Python -> a[39] in Julia
    xappa = (pdyn / 2.0)^a[39]
    STATE[].rh0 = a[40]; STATE[].g = a[41]  # a[39], a[40] in Python
    xappa3 = xappa^3
    xx, yy, zz = x * xappa, y * xappa, z * xappa
    sps, cps = sin(ps), cos(ps)
    x0 = A0_X0 / xappa; am = A0_A / xappa; s0 = A0_S0

    theta = (byimf == 0.0 && bzimf == 0.0) ? 0.0 : (t = atan(byimf, bzimf); t <= 0 ? t + 2π : t)
    sthetah = sin(theta / 2.0)^2
    factimf = a[24] + a[25] * sthetah  # a[23], a[24] in Python
    oimfy, oimfz = byimf * factimf, bzimf * factimf

    r = sqrt(x^2 + y^2 + z^2); xss, zss = x, z
    for _ in 1:20
        rh = STATE[].rh0 + rh2 * (zss / r)^2
        sinpsas = sps / (1.0 + (r / rh)^3)^0.33333333
        cospsas = sqrt(1.0 - sinpsas^2)
        xss_new = x * cospsas - z * sinpsas
        zss_new = x * sinpsas + z * cospsas
        abs(xss_new - xss) + abs(zss_new - zss) < 1.0e-6 && break
        xss, zss = xss_new, zss_new
    end

    rho2 = y^2 + zss^2; asq = am^2
    xmxm = max(am + xss - x0, 0.0); axx0 = xmxm^2; aro = asq + rho2
    sigma = sqrt((aro + axx0 + sqrt((aro + axx0)^2 - 4.0 * asq * axx0)) / (2.0 * asq))
    dsig = 0.003

    if sigma < (s0 + dsig)
        bxcf, bycf, bzcf = shlcar3x3(xx, yy, zz, ps)
        bxcf *= xappa3; bycf *= xappa3; bzcf *= xappa3

        # Tail: a[25], a[26], a[27], a[28] in Python -> a[26], a[27], a[28], a[29] in Julia
        STATE[].dxshift1 = a[26] + a[27] * g2; STATE[].dxshift2 = 0.0
        STATE[].d = a[28]; STATE[].deltady = a[29]
        bxt1, byt1, bzt1, bxt2, byt2, bzt2 = deformed(ps, xx, yy, zz)

        # Birk: a[34], a[35], a[36], a[37] in Python -> a[35], a[36], a[37], a[38] in Julia
        STATE[].xkappa1 = a[35] + a[36] * g2; STATE[].xkappa2 = a[37] + a[38] * g2
        bxr11, byr11, bzr11, bxr12, byr12, bzr12, bxr21, byr21, bzr21, bxr22, byr22, bzr22 = birk_tot(ps, xx, yy, zz)

        # RC: a[33], a[29], a[30], a[31], a[32] in Python -> a[34], a[30], a[31], a[32], a[33] in Julia
        STATE[].phi = 0.5π * tanh(abs(dst) / a[34])
        znam = max(abs(dst), 20.0)
        STATE[].sc_sy = a[30] * (20.0 / znam)^a[31] * xappa
        STATE[].sc_pr = a[32] * (20.0 / znam)^a[33] * xappa
        bxsrc, bysrc, bzsrc, bxprc, byprc, bzprc = full_rc(ps, xx, yy, zz)

        # Amplitudes: a[41], a[42] in Python -> a[42], a[43] in Julia
        dlp1 = (pdyn / 2.0)^a[42]; dlp2 = (pdyn / 2.0)^a[43]
        # a[1]-a[8] in Python -> a[2]-a[9] in Julia
        tamp1 = a[2] + a[3] * dlp1 + a[4] * g1 + a[5] * dst
        tamp2 = a[6] + a[7] * dlp2 + a[8] * g1 + a[9] * dst
        # a[9]-a[14] in Python -> a[10]-a[15] in Julia
        a_src = a[10] + a[11] * dst + a[12] * sqrt(pdyn)
        a_prc = a[13] + a[14] * dst + a[15] * sqrt(pdyn)
        # a[15]-a[22] in Python -> a[16]-a[23] in Julia
        a_r11 = a[16] + a[17] * g2; a_r12 = a[18] + a[19] * g2
        a_r21 = a[20] + a[21] * g2; a_r22 = a[22] + a[23] * g2

        # a[0] in Python -> a[1] in Julia, a[23], a[24] in Python -> a[24], a[25] in Julia
        bbx = a[1] * bxcf + tamp1 * bxt1 + tamp2 * bxt2 + a_src * bxsrc + a_prc * bxprc + a_r11 * bxr11 + a_r12 * bxr12 + a_r21 * bxr21 + a_r22 * bxr22 + a[24] * 0.0 + a[25] * 0.0 * sthetah
        bby = a[1] * bycf + tamp1 * byt1 + tamp2 * byt2 + a_src * bysrc + a_prc * byprc + a_r11 * byr11 + a_r12 * byr12 + a_r21 * byr21 + a_r22 * byr22 + a[24] * byimf + a[25] * byimf * sthetah
        bbz = a[1] * bzcf + tamp1 * bzt1 + tamp2 * bzt2 + a_src * bzsrc + a_prc * bzprc + a_r11 * bzr11 + a_r12 * bzr12 + a_r21 * bzr21 + a_r22 * bzr22 + a[24] * bzimf + a[25] * bzimf * sthetah

        if sigma < (s0 - dsig)
            return bbx, bby, bbz
        else
            fint = 0.5 * (1.0 - (sigma - s0) / dsig); fext = 1.0 - fint
            qx, qy, qz = t01_dipole(ps, x, y, z)
            return (bbx + qx) * fint + 0.0 * fext - qx, (bby + qy) * fint + oimfy * fext - qy, (bbz + qz) * fint + oimfz * fext - qz
        end
    else
        qx, qy, qz = t01_dipole(ps, x, y, z)
        return -qx, oimfy - qy, oimfz - qz
    end
end

function t01_dipole(ps, x, y, z)
    return dipole(ps, x, y, z; q0 = 30115.0)
end

function shlcar3x3(x, y, z, ps)
    a = SHLCAR_A
    p1, p2, p3 = a[37], a[38], a[39]
    r1, r2, r3 = a[40], a[41], a[42]
    q1, q2, q3 = a[43], a[44], a[45]
    s1, s2, s3 = a[46], a[47], a[48]
    t1, t2 = a[49], a[50]

    cps, sps = cos(ps), sin(ps)
    s2ps = 2.0 * cps
    st1, ct1 = sin(ps * t1), cos(ps * t1)
    st2, ct2 = sin(ps * t2), cos(ps * t2)
    x1, z1 = x * ct1 - z * st1, x * st1 + z * ct1
    x2, z2 = x * ct2 - z * st2, x * st2 + z * ct2

    # First sum (perpendicular symmetry) - 9 terms for p1,p2,p3 x r1,r2,r3
    p_arr = (p1, p2, p3)
    r_arr = (r1, r2, r3)

    hx_arr = zeros(9); hy_arr = zeros(9); hz_arr = zeros(9)
    idx = 0
    for pi in p_arr
        cyp, syp = cos(y / pi), sin(y / pi)
        for (k, rk) in enumerate(r_arr)
            idx += 1
            czr, szr = cos(z1 / rk), sin(z1 / rk)
            sqpr = sqrt(1.0 / pi^2 + 1.0 / rk^2)
            expr = exp(sqpr * x1)
            if k < 3
                fx = -sqpr * expr * cyp * szr
                hy_arr[idx] = expr / pi * syp * szr
                fz = -expr * cyp / rk * czr
            else
                # Special formula for k=3 (r3)
                fx = -expr * cyp * (sqpr * z1 * czr + szr / rk * (x1 + 1.0 / sqpr))
                hy_arr[idx] = expr / pi * syp * (z1 * czr + x1 / rk * szr / sqpr)
                fz = -expr * cyp * (czr * (1.0 + x1 / rk^2 / sqpr) - z1 / rk * szr)
            end
            hx_arr[idx] = fx * ct1 + fz * st1
            hz_arr[idx] = -fx * st1 + fz * ct1
        end
    end

    # Combine with coefficients
    bx = (a[1] + a[2] * cps) * hx_arr[1] + (a[3] + a[4] * cps) * hx_arr[2] + (a[5] + a[6] * cps) * hx_arr[3] +
        (a[7] + a[8] * cps) * hx_arr[4] + (a[9] + a[10] * cps) * hx_arr[5] + (a[11] + a[12] * cps) * hx_arr[6] +
        (a[13] + a[14] * cps) * hx_arr[7] + (a[15] + a[16] * cps) * hx_arr[8] + (a[17] + a[18] * cps) * hx_arr[9]
    by = (a[1] + a[2] * cps) * hy_arr[1] + (a[3] + a[4] * cps) * hy_arr[2] + (a[5] + a[6] * cps) * hy_arr[3] +
        (a[7] + a[8] * cps) * hy_arr[4] + (a[9] + a[10] * cps) * hy_arr[5] + (a[11] + a[12] * cps) * hy_arr[6] +
        (a[13] + a[14] * cps) * hy_arr[7] + (a[15] + a[16] * cps) * hy_arr[8] + (a[17] + a[18] * cps) * hy_arr[9]
    bz = (a[1] + a[2] * cps) * hz_arr[1] + (a[3] + a[4] * cps) * hz_arr[2] + (a[5] + a[6] * cps) * hz_arr[3] +
        (a[7] + a[8] * cps) * hz_arr[4] + (a[9] + a[10] * cps) * hz_arr[5] + (a[11] + a[12] * cps) * hz_arr[6] +
        (a[13] + a[14] * cps) * hz_arr[7] + (a[15] + a[16] * cps) * hz_arr[8] + (a[17] + a[18] * cps) * hz_arr[9]

    # Second sum (parallel symmetry) - 9 terms for q1,q2,q3 x s1,s2,s3
    q_arr = (q1, q2, q3)
    s_arr = (s1, s2, s3)

    idx = 0
    for qi in q_arr
        cyq, syq = cos(y / qi), sin(y / qi)
        for sk in s_arr
            idx += 1
            czs, szs = cos(z2 / sk), sin(z2 / sk)
            sqqs = sqrt(1.0 / qi^2 + 1.0 / sk^2)
            exqs = exp(sqqs * x2)
            fx = -sqqs * exqs * cyq * czs * sps
            hy = exqs / qi * syq * czs * sps
            fz = exqs * cyq / sk * szs * sps
            hx = fx * ct2 + fz * st2
            hz = -fx * st2 + fz * ct2
            coef = a[18 + 2 * idx - 1] + a[18 + 2 * idx] * s2ps
            bx += coef * hx
            by += coef * hy
            bz += coef * hz
        end
    end

    return bx, by, bz
end

include("t01_tail.jl")
include("t01_rc.jl")
include("t01_birk.jl")
