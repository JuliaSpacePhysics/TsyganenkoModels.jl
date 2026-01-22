# T01 Core Functions
using ..TsyganenkoModels: dipole
using LinearAlgebra: dot

mutable struct T01State
    dxshift1::Float64
    dxshift2::Float64
    d::Float64
    deltady::Float64
    sc_sy::Float64
    sc_pr::Float64
    g::Float64
end
const STATE = Ref(T01State(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0))

function t01_extall(a, pdyn, dst, byimf, bzimf, g1, g2, ps, x, y, z)
    rh2 = -5.2
    xappa = (pdyn / 2.0)^a[39]
    rh0 = a[40]
    STATE[].g = a[41]
    xappa3 = xappa^3
    xx, yy, zz = x * xappa, y * xappa, z * xappa
    sps = sin(ps)
    x0 = A0_X0 / xappa; am = A0_A / xappa; s0 = A0_S0
    bimf = (0, byimf, bzimf)

    theta = (byimf == 0.0 && bzimf == 0.0) ? 0.0 : (t = atan(byimf, bzimf); t <= 0 ? t + 2π : t)
    sthetah = sin(theta / 2.0)^2
    factimf = a[24] + a[25] * sthetah
    oimfy, oimfz = byimf * factimf, bzimf * factimf
    oimf = (0, oimfy, oimfz)

    r = sqrt(x^2 + y^2 + z^2); xss, zss = x, z
    for _ in 1:20
        rh = rh0 + rh2 * (zss / r)^2
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

    return if sigma < (s0 + dsig)
        bcf = shlcar3x3(xx, yy, zz, ps)

        # Tail
        STATE[].dxshift1 = a[26] + a[27] * g2; STATE[].dxshift2 = 0.0
        STATE[].d = a[28]; STATE[].deltady = a[29]
        bt1, bt2 = deformed(ps, xx, yy, zz, rh0)

        # Birk
        xkappa1 = a[35] + a[36] * g2; xkappa2 = a[37] + a[38] * g2
        br11, br12, br21, br22 = birk_tot(ps, xx, yy, zz, xkappa1, xkappa2)

        # RC
        phi = 0.5π * tanh(abs(dst) / a[34])
        znam = max(abs(dst), 20.0)
        STATE[].sc_sy = a[30] * (20.0 / znam)^a[31] * xappa
        STATE[].sc_pr = a[32] * (20.0 / znam)^a[33] * xappa
        bsrc, bprc = full_rc(ps, xx, yy, zz, phi)

        # Amplitudes
        dlp1 = (pdyn / 2.0)^a[42]; dlp2 = (pdyn / 2.0)^a[43]
        tamp1 = a[2] + a[3] * dlp1 + a[4] * g1 + a[5] * dst
        tamp2 = a[6] + a[7] * dlp2 + a[8] * g1 + a[9] * dst
        a_src = a[10] + a[11] * dst + a[12] * sqrt(pdyn)
        a_prc = a[13] + a[14] * dst + a[15] * sqrt(pdyn)
        a_r11 = a[16] + a[17] * g2; a_r12 = a[18] + a[19] * g2
        a_r21 = a[20] + a[21] * g2; a_r22 = a[22] + a[23] * g2

        b = @. a[1] * xappa3 * bcf + tamp1 * bt1 + tamp2 * bt2 + a_src * bsrc + a_prc * bprc + a_r11 * br11 + a_r12 * br12 + a_r21 * br21 + a_r22 * br22 + a[24] * bimf + a[25] * bimf * sthetah
        if sigma < (s0 - dsig)
            return b
        else
            fint = 0.5 * (1.0 - (sigma - s0) / dsig); fext = 1.0 - fint
            q = t01_dipole(ps, x, y, z)
            return @. (b + q) * fint + oimf * fext - q
        end
    else
        oimf .- t01_dipole(ps, x, y, z)
    end
end

function t01_dipole(ps, x, y, z)
    return dipole(ps, x, y, z; q0 = 30115.0)
end

@views function shlcar3x3(x, y, z, ps)
    a = SHLCAR_A
    p1, p2, p3 = a[37], a[38], a[39]
    r1, r2, r3 = a[40], a[41], a[42]
    q1, q2, q3 = a[43], a[44], a[45]
    s1, s2, s3 = a[46], a[47], a[48]
    t1, t2 = a[49], a[50]

    sps, cps = sincos(ps)
    s2ps = 2.0 * cps
    st1, ct1 = sincos(ps * t1)
    st2, ct2 = sincos(ps * t2)
    x1, z1 = x * ct1 - z * st1, x * st1 + z * ct1
    x2, z2 = x * ct2 - z * st2, x * st2 + z * ct2

    # First sum (perpendicular symmetry) - 9 terms for p1,p2,p3 x r1,r2,r3
    p_arr = (p1, p2, p3)
    r_arr = (r1, r2, r3)

    hx_arr = zeros(9); hy_arr = zeros(9); hz_arr = zeros(9)
    idx = 0
    for pi in p_arr
        syp, cyp = sincos(y / pi)
        for (k, rk) in enumerate(r_arr)
            idx += 1
            szr, czr = sincos(z1 / rk)
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

    # Combine with coefficients using broadcasting
    coef_perp = @. a[1:2:17] + a[2:2:18] * cps
    bx = dot(coef_perp, hx_arr)
    by = dot(coef_perp, hy_arr)
    bz = dot(coef_perp, hz_arr)

    # Second sum (parallel symmetry) - 9 terms for q1,q2,q3 x s1,s2,s3
    q_arr = (q1, q2, q3)
    s_arr = (s1, s2, s3)

    idx = 0
    for qi in q_arr
        syq, cyq = sincos(y / qi)
        for sk in s_arr
            idx += 1
            szs, czs = sincos(z2 / sk)
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
