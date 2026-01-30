function extall(pdyn, dst, byimf, bzimf, w1, w2, w3, w4, w5, w6, ps, x, y, z)
    a = TS04_A
    rh2 = -5.2
    dsig = 0.005  # T04 uses 0.005 vs T01's 0.003

    # T04 uses fixed rh0 and g values (different from T01)
    rh0 = 7.5   # tail hinging distance (T01 uses a[40] ≈ 9.0)
    g = 35.0    # tail warping parameter (T01 uses a[41] ≈ 28.2)

    xappa = (pdyn / 2.0)^a[23]
    xappa3 = xappa^3
    xx, yy, zz = x * xappa, y * xappa, z * xappa

    sps = sin(ps)
    x0 = TS04_A0_X0 / xappa
    am = TS04_A0_A / xappa
    s0 = TS04_A0_S0

    # IMF penetration factor - T04 uses single value (T01 uses theta-dependent)
    factimf = a[20]
    oimfy, oimfz = byimf * factimf, bzimf * factimf
    oimf = (0.0, oimfy, oimfz)

    # Iterative search for sigma (magnetopause distance)
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

    rho2 = y^2 + zss^2
    asq = am^2
    xmxm = max(am + xss - x0, 0.0)
    axx0 = xmxm^2
    aro = asq + rho2
    sigma = sqrt((aro + axx0 + sqrt((aro + axx0)^2 - 4.0 * asq * axx0)) / (2.0 * asq))

    return if sigma < (s0 + dsig)
        # Dipole shielding field
        bcf = shlcar3x3(xx, yy, zz, ps)

        # Tail field parameters - T04 uses dst-dependent dxshift and w1-dependent d
        dstt = min(dst, -20.0)
        znam_tail = abs(dstt)^0.37
        dxshift1 = a[24] - a[25] / znam_tail  # a[23]-a[24]/znam in Python
        dxshift2 = a[26] - a[27] / znam_tail  # a[25]-a[26]/znam in Python
        d = a[36] * exp(-w1 / a[37]) + a[69]  # a[35]*exp(-w1/a[36])+a[68] in Python
        deltady = 4.7  # T04 uses fixed value (T01 uses a[29])

        state = (; dxshift1, dxshift2, d, deltady, g)
        bt1, bt2 = deformed(ps, xx, yy, zz, rh0, state)

        # Birkeland field parameters - T04 uses dst-dependent xkappa
        znam_birk = max(abs(dst), 20.0)
        xkappa1 = a[32] * (znam_birk / 20.0)^a[33]
        xkappa2 = a[34] * (znam_birk / 20.0)^a[35]
        br11, _, br21, _ = birk_tot(ps, xx, yy, zz, xkappa1, xkappa2)

        # Ring current parameters - T04 uses fixed phi from array
        phi = a[38]
        znam_rc = max(abs(dst), 20.0)
        sc_sy = a[28] * (20.0 / znam_rc)^a[29] * xappa
        sc_pr = a[30] * (20.0 / znam_rc)^a[31] * xappa
        bsrc, bprc = full_rc(ps, xx, yy, zz, phi, sc_sy, sc_pr)

        # IMF components inside magnetosphere
        bimf = (0.0, byimf, bzimf)

        # Amplitude formulas - T04 uses w-indices with saturation terms
        dlp1 = (pdyn / 2.0)^a[21]
        dlp2 = (pdyn / 2.0)^a[22]

        # Saturation function: a*w/sqrt(w^2+a^2) approaches 1 as w→∞
        tamp1 = a[2] + a[3] * dlp1 + a[4] * a[39] * w1 / sqrt(w1^2 + a[39]^2) + a[5] * dst
        tamp2 = a[6] + a[7] * dlp2 + a[8] * a[40] * w2 / sqrt(w2^2 + a[40]^2) + a[9] * dst
        a_src = a[10] + a[11] * a[41] * w3 / sqrt(w3^2 + a[41]^2) + a[12] * dst
        a_prc = a[13] + a[14] * a[42] * w4 / sqrt(w4^2 + a[42]^2) + a[15] * dst
        a_r11 = a[16] + a[17] * a[43] * w5 / sqrt(w5^2 + a[43]^2)
        a_r21 = a[18] + a[19] * a[44] * w6 / sqrt(w6^2 + a[44]^2)

        # T04 only uses r11 and r21 Birkeland terms (not r12 and r22 like T01)
        b = @. a[1] * xappa3 * bcf + tamp1 * bt1 + tamp2 * bt2 +
            a_src * bsrc + a_prc * bprc +
            a_r11 * br11 + a_r21 * br21 +
            a[20] * bimf

        if sigma < (s0 - dsig)
            return b
        else
            # Interpolation region near magnetopause
            fint = 0.5 * (1.0 - (sigma - s0) / dsig)
            fext = 1.0 - fint
            q = dipole(ps, x, y, z; q0 = 30115.0)
            return @. (b + q) * fint + oimf * fext - q
        end
    else
        # Outside magnetosphere
        oimf .- dipole(ps, x, y, z; q0 = 30115.0)
    end
end