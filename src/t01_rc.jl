# T01 Ring Current Functions

function full_rc(ps, x, y, z)
    hxsrc, hysrc, hzsrc, hxprc, hyprc, hzprc = src_prc(ps, x, y, z)
    x_sc = STATE[].sc_sy - 1.0
    fsx, fsy, fsz = rc_shield(C_SY, ps, x_sc, x, y, z)
    x_sc = STATE[].sc_pr - 1.0
    fpx, fpy, fpz = rc_shield(C_PR, ps, x_sc, x, y, z)
    return hxsrc + fsx, hysrc + fsy, hzsrc + fsz, hxprc + fpx, hyprc + fpy, hzprc + fpz
end

function src_prc(ps, x, y, z)
    cps, sps = cos(ps), sin(ps)
    xt, zt = x * cps - z * sps, z * cps + x * sps
    xts, yts, zts = xt / STATE[].sc_sy, y / STATE[].sc_sy, zt / STATE[].sc_sy
    xta, yta, zta = xt / STATE[].sc_pr, y / STATE[].sc_pr, zt / STATE[].sc_pr
    bxs, bys, bzs = rc_symm(xts, yts, zts)
    bxa_s, bya_s, bza_s = prc_symm(xta, yta, zta)
    cp, sp = cos(STATE[].phi), sin(STATE[].phi)
    xr, yr = xta * cp - yta * sp, xta * sp + yta * cp
    bxa_qr, bya_qr, bza_q = prc_quad(xr, yr, zta)
    bxa_q, bya_q = bxa_qr * cp + bya_qr * sp, -bxa_qr * sp + bya_qr * cp
    bxp, byp, bzp = bxa_s + bxa_q, bya_s + bya_q, bza_s + bza_q
    bxsrc, bysrc, bzsrc = bxs * cps + bzs * sps, bys, bzs * cps - bxs * sps
    bxprc, byprc, bzprc = bxp * cps + bzp * sps, byp, bzp * cps - bxp * sps
    return bxsrc, bysrc, bzsrc, bxprc, byprc, bzprc
end

function rc_symm(x, y, z)
    ds, dc, d, drd = 1.0e-2, 0.99994999875, 1.0e-4, 5.0e3
    rho2 = x^2 + y^2; r2 = rho2 + z^2; r = sqrt(r2)
    rp, rm = r + d, r - d; sint, cost = sqrt(rho2) / r, z / r
    if sint < ds
        a = ap_rc(r, ds, dc) / ds
        dardr = (rp * ap_rc(rp, ds, dc) - rm * ap_rc(rm, ds, dc)) * drd
        fxy = z * (2.0 * a - dardr) / (r * r2)
        return fxy * x, fxy * y, (2.0 * a * cost^2 + dardr * sint^2) / r
    else
        theta = atan(sint, cost); tp, tm = theta + d, theta - d
        sintp, costp, sintm, costm = sin(tp), cos(tp), sin(tm), cos(tm)
        br = (sintp * ap_rc(r, sintp, costp) - sintm * ap_rc(r, sintm, costm)) / (r * sint) * drd
        bt = (rm * ap_rc(rm, sint, cost) - rp * ap_rc(rp, sint, cost)) / r * drd
        fxy = (br + bt * cost / sint) / r
        return fxy * x, fxy * y, br * cost - bt * sint
    end
end

function ap_rc(r, sint, cost)
    a1, a2 = -456.5289941, 375.9055332
    rrc1, dd1, rrc2, dd2 = 4.27468495, 2.439528329, 3.367557287, 3.146382545
    p1, r1, dr1, dla1 = -0.2291904607, 3.74606474, 1.508802177, 0.5873525737
    p2, r2, dr2, dla2 = 0.1556236119, 4.993638842, 3.324180497, 0.4368407663
    p3, r3, dr3 = 0.1855957207, 2.969226745, 2.243367377
    prox = false; sint1, cost1 = sint, cost
    if sint1 < 1.0e-2
        sint1, cost1, prox = 1.0e-2, 0.99994999875, true
    end
    alpha, gamma = sint1^2 / r, cost1 / r^2
    arg1 = -((r - r1) / dr1)^2 - (cost1 / dla1)^2
    arg2 = -((r - r2) / dr2)^2 - (cost1 / dla2)^2
    arg3 = -((r - r3) / dr3)^2
    dexp1 = arg1 < -500.0 ? 0.0 : exp(arg1)
    dexp2 = arg2 < -500.0 ? 0.0 : exp(arg2)
    dexp3 = arg3 < -500.0 ? 0.0 : exp(arg3)
    alpha_s = alpha * (1.0 + p1 * dexp1 + p2 * dexp2 + p3 * dexp3)
    gammas2 = gamma^2; alsqh = alpha_s^2 / 2.0
    f = 64.0 / 27.0 * gammas2 + alsqh^2
    q = (sqrt(f) + alsqh)^(1.0 / 3.0)
    c = max(q - 4.0 * gammas2^(1.0 / 3.0) / (3.0 * q), 0.0)
    g = sqrt(c^2 + 4.0 * gammas2^(1.0 / 3.0))
    rs = 4.0 / ((sqrt(2.0 * g - c) + sqrt(c)) * (g + c))
    costs, sints = gamma * rs^2, sqrt(1.0 - (gamma * rs^2)^2)
    rhos, zs = rs * sints, rs * costs
    ap = a1 * elliptic_aphi(rrc1, rhos, zs, dd1) + a2 * elliptic_aphi(rrc2, rhos, zs, dd2)
    return prox ? ap * sint / sint1 : ap
end

function elliptic_aphi(rrc, rhos, zs, dd)
    p = (rrc + rhos)^2 + zs^2 + dd^2; xk2 = 4.0 * rrc * rhos / p
    xk = sqrt(xk2); xkrho12 = xk * sqrt(rhos); xk2s = 1.0 - xk2; dl = log(1.0 / xk2s)
    elk = 1.38629436112 + xk2s * (0.09666344259 + xk2s * (0.03590092383 + xk2s * (0.03742563713 + xk2s * 0.01451196212))) + dl * (0.5 + xk2s * (0.12498593597 + xk2s * (0.06880248576 + xk2s * (0.03328355346 + xk2s * 0.00441787012))))
    ele = 1.0 + xk2s * (0.44325141463 + xk2s * (0.0626060122 + xk2s * (0.04757383546 + xk2s * 0.01736506451))) + dl * xk2s * (0.2499836831 + xk2s * (0.09200180037 + xk2s * (0.04069697526 + xk2s * 0.00526449639)))
    return ((1.0 - xk2 * 0.5) * elk - ele) / xkrho12
end

function prc_symm(x, y, z)
    ds, dc, d, drd = 1.0e-2, 0.99994999875, 1.0e-4, 5.0e3
    rho2 = x^2 + y^2; r2 = rho2 + z^2; r = sqrt(r2)
    rp, rm = r + d, r - d; sint, cost = sqrt(rho2) / r, z / r
    if sint < ds
        a = apprc(r, ds, dc) / ds
        dardr = (rp * apprc(rp, ds, dc) - rm * apprc(rm, ds, dc)) * drd
        fxy = z * (2.0 * a - dardr) / (r * r2)
        return fxy * x, fxy * y, (2.0 * a * cost^2 + dardr * sint^2) / r
    else
        theta = atan(sint, cost); tp, tm = theta + d, theta - d
        sintp, costp, sintm, costm = sin(tp), cos(tp), sin(tm), cos(tm)
        br = (sintp * apprc(r, sintp, costp) - sintm * apprc(r, sintm, costm)) / (r * sint) * drd
        bt = (rm * apprc(rm, sint, cost) - rp * apprc(rp, sint, cost)) / r * drd
        fxy = (br + bt * cost / sint) / r
        return fxy * x, fxy * y, br * cost - bt * sint
    end
end

function apprc(r, sint, cost)
    a1, a2 = -80.11202281, 12.58246758
    rrc1, dd1, rrc2, dd2 = 6.560486035, 1.930711037, 3.827208119, 0.7789990504
    p1, alpha1, dal1, beta1, dg1 = 0.3058309043, 0.1817139853, 0.1257532909, 3.422509402, 0.04742939676
    p2, alpha2, dal2, beta2, dg2, beta3 = -4.800458958, -0.02845643596, 0.2188114228, 2.545944574, 0.00813272793, 0.35868244
    p3, alpha3, dal3, beta4, dg3, beta5 = 103.1601001, -0.00764731187, 0.1046487459, 2.958863546, 0.01172314188, 0.4382872938
    q0, q1, alpha4, dal4, dg4 = 0.0113490815, 14.51339943, 0.2647095287, 0.07091230197, 0.01512963586
    q2, alpha5, dal5, dg5, beta6, beta7 = 6.861329631, 0.1677400816, 0.04433648846, 0.05553741389, 0.7665599464, 0.7277854652
    prox = false; sint1, cost1 = sint, cost
    if sint1 < 1.0e-2
        sint1, cost1, prox = 1.0e-2, 0.99994999875, true
    end
    alpha, gamma = sint1^2 / r, cost1 / r^2
    arg1 = -(gamma / dg1)^2
    arg2 = -((alpha - alpha4) / dal4)^2 - (gamma / dg4)^2
    dexp1 = arg1 < -500.0 ? 0.0 : exp(arg1)
    dexp2 = arg2 < -500.0 ? 0.0 : exp(arg2)
    alpha_s = alpha * (1.0 + p1 / (1.0 + ((alpha - alpha1) / dal1)^2)^beta1 * dexp1 + p2 * (alpha - alpha2) / (1.0 + ((alpha - alpha2) / dal2)^2)^beta2 / (1.0 + (gamma / dg2)^2)^beta3 + p3 * (alpha - alpha3)^2 / (1.0 + ((alpha - alpha3) / dal3)^2)^beta4 / (1.0 + (gamma / dg3)^2)^beta5)
    gamma_s = gamma * (1.0 + q0 + q1 * (alpha - alpha4) * dexp2 + q2 * (alpha - alpha5) / (1.0 + ((alpha - alpha5) / dal5)^2)^beta6 / (1.0 + (gamma / dg5)^2)^beta7)
    gammas2 = gamma_s^2; alsqh = alpha_s^2 / 2.0
    f = 64.0 / 27.0 * gammas2 + alsqh^2
    q_val = (sqrt(f) + alsqh)^(1.0 / 3.0)
    c = max(q_val - 4.0 * gammas2^(1.0 / 3.0) / (3.0 * q_val), 0.0)
    g = sqrt(c^2 + 4.0 * gammas2^(1.0 / 3.0))
    rs = 4.0 / ((sqrt(2.0 * g - c) + sqrt(c)) * (g + c))
    costs, sints = gamma_s * rs^2, sqrt(1.0 - (gamma_s * rs^2)^2)
    rhos, zs = rs * sints, rs * costs
    ap = a1 * elliptic_aphi(rrc1, rhos, zs, dd1) + a2 * elliptic_aphi(rrc2, rhos, zs, dd2)
    return prox ? ap * sint / sint1 : ap
end

function prc_quad(x, y, z)
    d, dd, ds, dc = 1.0e-4, 2.0e-4, 1.0e-2, 0.99994999875
    rho2 = x^2 + y^2; r = sqrt(rho2 + z^2); rho = sqrt(rho2)
    sint, cost = rho / r, z / r; rp, rm = r + d, r - d
    if sint > ds
        cphi, sphi = x / rho, y / rho
        br, bt = br_prc_q(r, sint, cost), bt_prc_q(r, sint, cost)
        dbrr = (br_prc_q(rp, sint, cost) - br_prc_q(rm, sint, cost)) / dd
        theta = atan(sint, cost); tp, tm = theta + d, theta - d
        sintp, costp, sintm, costm = sin(tp), cos(tp), sin(tm), cos(tm)
        dbtt = (bt_prc_q(r, sintp, costp) - bt_prc_q(r, sintm, costm)) / dd
        bx = sint * (br + (br + r * dbrr + dbtt) * sphi^2) + cost * bt
        by = -sint * sphi * cphi * (br + r * dbrr + dbtt)
        bz = (br * cost - bt * sint) * cphi
    else
        st, ct = ds, z < 0.0 ? -dc : dc
        theta = atan(st, ct); tp, tm = theta + d, theta - d
        sintp, costp, sintm, costm = sin(tp), cos(tp), sin(tm), cos(tm)
        br, bt = br_prc_q(r, st, ct), bt_prc_q(r, st, ct)
        dbrr = (br_prc_q(rp, st, ct) - br_prc_q(rm, st, ct)) / dd
        dbtt = (bt_prc_q(r, sintp, costp) - bt_prc_q(r, sintm, costm)) / dd
        fcxy = r * dbrr + dbtt
        bx = (br * (x^2 + 2.0 * y^2) + fcxy * y^2) / (r * st)^2 + bt * cost
        by = -(br + fcxy) * x * y / (r * st)^2
        bz = (br * cost / st - bt) * x / r
    end
    return bx, by, bz
end

function ffs(a, a0, da)
    sq1 = sqrt((a + a0)^2 + da^2); sq2 = sqrt((a - a0)^2 + da^2)
    fa = 2.0 / (sq1 + sq2); f = fa * a
    fs = 0.5 * (sq1 + sq2) / (sq1 * sq2) * (1.0 - f^2)
    return f, fa, fs
end

function br_prc_q(r, sint, cost)
    a = (-21.2666329, 32.24527521, -6.062894078, 7.515660734, 233.7341288, -227.1195714, 8.483233889, 16.80642754, -24.63534184, 9.067120578, -1.052686913, -12.08384538, 18.61969572, -12.71686069, 47017.35679, -50646.71204, 7746.058231, 1.531069371)
    xk1, al1, dal1, b1, be1 = 2.318824273, 0.1417519429, 0.638801311e-2, 5.303934488, 4.213397467
    xk2, al2, dal2, b2, be2 = 0.7955534018, 0.1401142771, 0.2306094179e-1, 3.462235072, 2.56874301
    xk3, xk4, al3, dal3, b3, be3 = 3.477425908, 1.92215511, 0.1485233485, 0.2319676273e-1, 7.830223587, 8.492933868
    al4, dal4, dg1 = 0.1295221828, 0.01753008801, 0.01125504083
    al5, dal5, dg2 = 0.1811846095, 0.04841237481, 0.01981805097
    c1, c2, c3 = 6.557801891, 6.348576071, 5.744436687
    al6, dal6, drm = 0.2265212965, 0.1301957209, 0.5654023158
    sint2, cost2, sc = sint^2, cost^2, sint * cost
    alpha, gamma = sint2 / r, cost / r^2
    f, fa, fs = ffs(alpha, al1, dal1); d1 = sc * f^xk1 / ((r / b1)^be1 + 1.0); d2 = d1 * cost2
    f, fa, fs = ffs(alpha, al2, dal2); d3 = sc * fs^xk2 / ((r / b2)^be2 + 1.0); d4 = d3 * cost2
    f, fa, fs = ffs(alpha, al3, dal3); d5 = sc * (alpha^xk3) * (fs^xk4) / ((r / b3)^be3 + 1.0); d6 = d5 * cost2
    arga = ((alpha - al4) / dal4)^2 + 1.0; argg = 1.0 + (gamma / dg1)^2
    d7 = sc / arga / argg; d8 = d7 / arga; d9 = d8 / arga; d10 = d9 / arga
    arga = ((alpha - al5) / dal5)^2 + 1.0; argg = 1.0 + (gamma / dg2)^2
    d11 = sc / arga / argg; d12 = d11 / arga; d13 = d12 / arga; d14 = d13 / arga
    d15 = sc / (r^4 + c1^4); d16 = sc / (r^4 + c2^4) * cost2; d17 = sc / (r^4 + c3^4) * cost2^2
    f, fa, fs = ffs(alpha, al6, dal6); d18 = sc * fs / (1.0 + ((r - 1.2) / drm)^2)
    return sum(a .* (d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12, d13, d14, d15, d16, d17, d18))
end

function bt_prc_q(r, sint, cost)
    a = (12.74640393, -7.516393516, -5.476233865, 3.212704645, -59.10926169, 46.62198189, -0.01644280062, 0.1234229112, -0.08579198697, 0.01321366966, 0.8970494003, 9.136186247, -38.19301215, 21.73775846, -410.0783424, -69.9083269, -848.854344)
    xk1, al1, dal1, b1, be1 = 1.243288286, 0.207172136, 0.05030555417, 7.471332374, 3.180533613
    xk2, al2, dal2, be2 = 1.376743507, 0.1568504222, 0.02092910682, 1.985148197
    xk3, xk4, al3, dal3, b3, be3 = 0.315713994, 1.056309517, 0.1701395257, 0.101987007, 6.293740981, 5.671824276
    al4, dal4, dg1 = 0.1280772299, 0.02189060799, 0.0104069608
    al5, dal5, dg2 = 0.1648265607, 0.04701592613, 0.01526400086
    c1, c2, c3 = 12.88384229, 3.361775101, 23.44173897
    sint2, cost2 = sint^2, cost^2; alpha, gamma = sint2 / r, cost / r^2
    f, fa, fs = ffs(alpha, al1, dal1); d1 = f^xk1 / ((r / b1)^be1 + 1.0); d2 = d1 * cost2
    f, fa, fs = ffs(alpha, al2, dal2); d3 = fa^xk2 / r^be2; d4 = d3 * cost2
    f, fa, fs = ffs(alpha, al3, dal3); d5 = fs^xk3 * alpha^xk4 / ((r / b3)^be3 + 1.0); d6 = d5 * cost2
    f, fa, fs = ffs(gamma, 0.0, dg1); fcc = 1.0 + ((alpha - al4) / dal4)^2
    d7 = fs / fcc; d8 = d7 / fcc; d9 = d8 / fcc; d10 = d9 / fcc
    arg = 1.0 + ((alpha - al5) / dal5)^2
    d11 = 1.0 / arg / (1.0 + (gamma / dg2)^2); d12 = d11 / arg; d13 = d12 / arg; d14 = d13 / arg
    d15 = 1.0 / (r^4 + c1^2); d16 = cost2 / (r^4 + c2^2); d17 = cost2^2 / (r^4 + c3^2)
    return sum(a .* (d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12, d13, d14, d15, d16, d17))
end

function rc_shield(a, ps, x_sc, x, y, z)
    fac_sc = (x_sc + 1.0)^3; cps, sps = cos(ps), sin(ps); s3ps = 2.0 * cps
    pst1, pst2 = ps * a[85], ps * a[86]
    st1, ct1 = sin(pst1), cos(pst1); st2, ct2 = sin(pst2), cos(pst2)
    x1, z1 = x * ct1 - z * st1, x * st1 + z * ct1
    x2, z2 = x * ct2 - z * st2, x * st2 + z * ct2
    bx, by, bz = 0.0, 0.0, 0.0; l = 0
    for m in 1:2, i in 1:3
        p, q = a[72 + i], a[78 + i]
        cypi, sypi = cos(y / p), sin(y / p); cyqi, syqi = cos(y / q), sin(y / q)
        for k in 1:3
            r, s = a[75 + k], a[81 + k]
            szrk, czrk = sin(z1 / r), cos(z1 / r); czsk, szsk = cos(z2 / s), sin(z2 / s)
            sqpr, sqqs = sqrt(1.0 / p^2 + 1.0 / r^2), sqrt(1.0 / q^2 + 1.0 / s^2)
            epr, eqs = exp(x1 * sqpr), exp(x2 * sqqs)
            for n in 1:2, nn in 1:2
                if m == 1
                    fx = -sqpr * epr * cypi * szrk * fac_sc
                    fy = epr * sypi * szrk / p * fac_sc
                    fz = -epr * cypi * czrk / r * fac_sc
                    mult = n == 1 ? (nn == 1 ? 1.0 : x_sc) : (nn == 1 ? cps : cps * x_sc)
                else
                    fx = -sps * sqqs * eqs * cyqi * czsk * fac_sc
                    fy = sps / q * eqs * syqi * czsk * fac_sc
                    fz = sps / s * eqs * cyqi * szsk * fac_sc
                    mult = n == 1 ? (nn == 1 ? 1.0 : x_sc) : (nn == 1 ? s3ps : s3ps * x_sc)
                end
                hx, hy, hz = fx * mult, fy * mult, fz * mult
                ct, st = m == 1 ? (ct1, st1) : (ct2, st2)
                hxr, hzr = hx * ct + hz * st, -hx * st + hz * ct
                l += 1; bx += hxr * a[l]; by += hy * a[l]; bz += hzr * a[l]
            end
        end
    end
    return bx, by, bz
end
