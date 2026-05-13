function t01_extall(pdyn, dst, byimf, bzimf, g1, g2, ps, x, y, z)
    c = T01_C
    xappa = (pdyn / 2.0)^c.xappa_exp
    xappa3 = xappa^3
    xx, yy, zz = x * xappa, y * xappa, z * xappa
    sps = sin(ps)
    x0 = MP_X0 / xappa; am = MP_A / xappa; s0 = MP_S0
    bimf = (0.0, byimf, bzimf)

    theta = (byimf == 0.0 && bzimf == 0.0) ? 0.0 : (t = atan(byimf, bzimf); t <= 0 ? t + 2π : t)
    sthetah = sin(theta / 2.0)^2
    factimf = c.factimf.c + c.factimf.sth * sthetah
    oimf = (0.0, byimf * factimf, bzimf * factimf)

    xss, zss = _warped_xz(x, y, z, sps, c.rh0)
    rho2 = y^2 + zss^2
    sigma = _sigma(xss, x0, am, rho2)

    return _switch(sigma, s0, 0.003, ps, x, y, z, oimf; q0 = 30115.0) do
        _tail_amp(c, dlp, g1, dst) = c.c + c.dlp * dlp + c.g1 * g1 + c.dst * dst
        _birk_amp(c) = c.c + c.g2 * g2
        _rc_amp(c) = c.c + c.dst * dst + c.sqp * sqrt(pdyn)

        bcf = shlcar3x3(xx, yy, zz, ps)
        # Tail
        state = (;
            dxshift1 = c.dxshift.c + c.dxshift.g2 * g2,
            dxshift2 = 0.0,
            d = c.d, deltady = c.deltady, g = c.g,
        )
        bt1, bt2 = deformed(ps, xx, yy, zz, c.rh0, state)
        dlp1 = (pdyn / 2.0)^c.dlp1_exp
        dlp2 = (pdyn / 2.0)^c.dlp2_exp
        B_tail = _tail_amp(c.tamp1, dlp1, g1, dst) .* bt1 .+ _tail_amp(c.tamp2, dlp2, g1, dst) .* bt2

        # Birkeland
        xkappa1 = c.xkappa1.c + c.xkappa1.g2 * g2
        xkappa2 = c.xkappa2.c + c.xkappa2.g2 * g2
        br11, br12, br21, br22 = birk_tot(ps, xx, yy, zz, xkappa1, xkappa2)
        B_birk = _birk_amp(c.r11) .* br11 .+ _birk_amp(c.r12) .* br12 .+ _birk_amp(c.r21) .* br21 .+ _birk_amp(c.r22) .* br22

        # Ring current
        phi = 0.5π * tanh(abs(dst) / c.rc_phi_denom)
        znam = max(abs(dst), 20.0)
        sc_sy = c.sc_sy.c * (20.0 / znam)^c.sc_sy.exp * xappa
        sc_pr = c.sc_pr.c * (20.0 / znam)^c.sc_pr.exp * xappa
        bsrc, bprc = full_rc(ps, xx, yy, zz, phi, sc_sy, sc_pr)
        B_rc = _rc_amp(c.src) .* bsrc .+ _rc_amp(c.prc) .* bprc

        @. c.amp_cf * xappa3 * bcf + B_tail + B_rc + B_birk + factimf * bimf
    end
end
