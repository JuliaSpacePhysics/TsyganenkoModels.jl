# T01 Core Functions
using LinearAlgebra: dot

function extall(pdyn, dst, byimf, bzimf, g1, g2, ps, x, y, z)
    a = T01_A
    rh2 = -5.2
    xappa = (pdyn / 2.0)^a[39]
    rh0 = a[40]
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
        state = (;
            dxshift1 = a[26] + a[27] * g2, dxshift2 = 0.0, d = a[28], deltady = a[29], g = a[41],
        )
        bt1, bt2 = deformed(ps, xx, yy, zz, rh0, state)

        # Birk
        xkappa1 = a[35] + a[36] * g2; xkappa2 = a[37] + a[38] * g2
        br11, br12, br21, br22 = birk_tot(ps, xx, yy, zz, xkappa1, xkappa2)

        # RC
        phi = 0.5π * tanh(abs(dst) / a[34])
        znam = max(abs(dst), 20.0)
        sc_sy = a[30] * (20.0 / znam)^a[31] * xappa
        sc_pr = a[32] * (20.0 / znam)^a[33] * xappa
        bsrc, bprc = full_rc(ps, xx, yy, zz, phi, sc_sy, sc_pr)

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
