# Saturation function from TS04: bounded surrogate for the W-index drivers.
@inline _sat(w, a) = a * w / sqrt(w^2 + a^2)

function ts04_extall(pdyn, dst, byimf, bzimf, w1, w2, w3, w4, w5, w6, ps, x, y, z)
    c = TS04_C
    xappa = (pdyn / 2.0)^c.xappa_exp
    xappa3 = xappa^3
    xx, yy, zz = x * xappa, y * xappa, z * xappa

    sps = sin(ps)
    x0 = MP_X0 / xappa; am = MP_A / xappa; s0 = MP_S0

    # TS04 uses a single (theta-independent) IMF penetration factor.
    oimf = (0.0, byimf * c.factimf, bzimf * c.factimf)

    xss, zss = _warped_xz(x, y, z, sps, c.rh0)
    rho2 = y^2 + zss^2
    sigma = _sigma(xss, x0, am, rho2)

    return _switch(sigma, s0, 0.005, ps, x, y, z, oimf; q0 = 30115.0) do
        _amp(c, w, β) = c.c + c.w * _sat(w, β)
        _tamp(c, dlp, w, β) = c.c + c.dlp * dlp + c.w * _sat(w, β)

        ws = c.w_sat

        bcf = shlcar3x3(xx, yy, zz, ps)

        # Tail (dst-dependent dxshift, w1-dependent thickness)
        znam_tail = abs(min(dst, -20.0))^0.37
        state = (;
            dxshift1 = c.dxshift1.c - c.dxshift1.znam / znam_tail,
            dxshift2 = c.dxshift2.c - c.dxshift2.znam / znam_tail,
            d = c.d.a * exp(-w1 / c.d.b) + c.d.off,
            deltady = c.deltady,
            g = c.g,
        )
        bt1, bt2 = deformed(ps, xx, yy, zz, c.rh0, state)
        dlp1 = (pdyn / 2.0)^c.dlp1_exp
        dlp2 = (pdyn / 2.0)^c.dlp2_exp
        B_tail = _tamp(c.tamp1, dlp1, w1, ws[1]) .* bt1 .+ _tamp(c.tamp2, dlp2, w2, ws[2]) .* bt2

        # Birkeland (TS04 keeps only region-1 of each mode; r12, r22 unused)
        znam = max(abs(dst), 20.0)
        xkappa1 = c.xkappa1.c * (znam / 20.0)^c.xkappa1.exp
        xkappa2 = c.xkappa2.c * (znam / 20.0)^c.xkappa2.exp
        br11, _, br21, _ = birk_tot(ps, xx, yy, zz, xkappa1, xkappa2)
        B_birk = _amp(c.r11, w5, ws[5]) .* br11 .+ _amp(c.r21, w6, ws[6]) .* br21

        # Ring current (fixed phi from fit)
        sc_sy = c.sc_sy.c * (20.0 / znam)^c.sc_sy.exp * xappa
        sc_pr = c.sc_pr.c * (20.0 / znam)^c.sc_pr.exp * xappa
        bsrc, bprc = full_rc(ps, xx, yy, zz, c.phi, sc_sy, sc_pr)
        B_rc = _amp(c.src, w3, ws[3]) .* bsrc .+ _amp(c.prc, w4, ws[4]) .* bprc

        # Amplitudes driven by W indices through the saturation function.
        bimf = (0.0, byimf, bzimf)


        @. c.amp_cf * xappa3 * bcf + B_tail + B_rc + B_birk + c.factimf * bimf
    end
end
