# T01 Birkeland Current Functions

function birk_tot(ps, x, y, z)
    x_sc11 = STATE[].xkappa1 - 1.1; x_sc21 = STATE[].xkappa2 - 1.0  # region1: -1.1, region2: -1.0
    fx11, fy11, fz11 = birk_1n2(1, 1, ps, x, y, z, STATE[].xkappa1)
    hx11, hy11, hz11 = birk_shl(SH11, ps, x_sc11, x, y, z)
    fx12, fy12, fz12 = birk_1n2(1, 2, ps, x, y, z, STATE[].xkappa1)
    hx12, hy12, hz12 = birk_shl(SH12, ps, x_sc11, x, y, z)
    fx21, fy21, fz21 = birk_1n2(2, 1, ps, x, y, z, STATE[].xkappa2)
    hx21, hy21, hz21 = birk_shl(SH21, ps, x_sc21, x, y, z)
    fx22, fy22, fz22 = birk_1n2(2, 2, ps, x, y, z, STATE[].xkappa2)
    hx22, hy22, hz22 = birk_shl(SH22, ps, x_sc21, x, y, z)
    return fx11 + hx11, fy11 + hy11, fz11 + hz11, fx12 + hx12, fy12 + hy12, fz12 + hz12, fx21 + hx21, fy21 + hy21, fz21 + hz21, fx22 + hx22, fy22 + hy22, fz22 + hz22
end

function birk_1n2(numb, mode, ps, x, y, z, xkappa)
    # Select coefficient array based on region and mode
    a = numb == 1 ? (mode == 1 ? A11 : A12) : (mode == 1 ? A21 : A22)

    # Parameters from Python birk_1n2
    beta, rh, eps, b_param, rho_0 = 0.9, 10.0, 3.0, 0.5, 7.0
    dphi = numb == 1 ? 0.055 : 0.03
    dtheta = numb == 1 ? 0.06 : 0.09

    # Scale coordinates
    xsc, ysc, zsc = x * xkappa, y * xkappa, z * xkappa
    rho = sqrt(xsc^2 + zsc^2)
    rsc = sqrt(xsc^2 + ysc^2 + zsc^2)
    rho2 = rho_0^2

    # Cylindrical phi
    phi = (xsc == 0.0 && zsc == 0.0) ? 0.0 : atan(-zsc, xsc)
    sphic, cphic = sin(phi), cos(phi)

    # Deformation
    brack = dphi + b_param * rho2 / (rho2 + 1.0) * (rho^2 - 1.0) / (rho2 + rho^2)
    r1rh = (rsc - 1.0) / rh
    psias = beta * ps / (1.0 + r1rh^eps)^(1.0 / eps)

    phis = phi - brack * sin(phi) - psias
    dphisphi = 1.0 - brack * cos(phi)
    dphisrho = -2.0 * b_param * rho2 * rho / (rho2 + rho^2)^2 * sin(phi) + beta * ps * r1rh^(eps - 1) * rho / (rh * rsc * (1.0 + r1rh^eps)^(1.0 / eps + 1))
    dphisdy = beta * ps * r1rh^(eps - 1) * ysc / (rh * rsc * (1.0 + r1rh^eps)^(1.0 / eps + 1))

    sphics, cphics = sin(phis), cos(phis)
    xs = rho * cphics
    zs = -rho * sphics

    # Call twocones with deformed coordinates
    bxs, byas, bzs = twocones(a, xs, ysc, zs, mode, dtheta)

    # Transform back
    brhoas = bxs * cphics - bzs * sphics
    bphias = -bxs * sphics - bzs * cphics

    brho_s = brhoas * dphisphi * xkappa
    bphi_s = (bphias - rho * (byas * dphisdy + brhoas * dphisrho)) * xkappa
    by_s = byas * dphisphi * xkappa

    bx = brho_s * cphic - bphi_s * sphic
    by = by_s
    bz = -brho_s * sphic - bphi_s * cphic

    return bx, by, bz
end

function twocones(a, x, y, z, mode, dtheta)
    bxn, byn, bzn = one_cone(a, x, y, z, mode, dtheta)
    bxs, bys, bzs = one_cone(a, x, -y, -z, mode, dtheta)
    return bxn - bxs, byn + bys, bzn + bzs
end

function one_cone(a, x, y, z, mode, dtheta)
    dr, dt = 1.0e-6, 1.0e-6
    theta0 = a[31]  # a[30] in Python (0-indexed)

    rho2 = x^2 + y^2; rho = sqrt(rho2)
    r = sqrt(rho2 + z^2)
    r < 1.0e-8 && return 0.0, 0.0, 0.0

    theta = atan(rho, z)
    phi = atan(y, x)

    # Deformation
    rs = r_s(a, r, theta)
    thetas = theta_s(a, r, theta)

    # Field at deformed position
    btast, bfast = fialcos(rs, thetas, phi, mode, theta0, dtheta)

    # Derivatives for transformation
    drsdr = (r_s(a, r + dr, theta) - r_s(a, r - dr, theta)) / (2.0 * dr)
    drsdt = (r_s(a, r, theta + dt) - r_s(a, r, theta - dt)) / (2.0 * dt)
    dtsdr = (theta_s(a, r + dr, theta) - theta_s(a, r - dr, theta)) / (2.0 * dr)
    dtsdt = (theta_s(a, r, theta + dt) - theta_s(a, r, theta - dt)) / (2.0 * dt)

    stsst = sin(thetas) / max(sin(theta), 1.0e-10)
    rsr = rs / r

    # Transform field (Python formula)
    br = -rsr / r * stsst * btast * drsdt
    btheta = rsr * stsst * btast * drsdr
    bphi = rsr * bfast * (drsdr * dtsdt - drsdt * dtsdr)

    s = rho / r
    c = z / r
    sf = rho > 1.0e-10 ? y / rho : 0.0
    cf = rho > 1.0e-10 ? x / rho : 1.0

    be = br * s + btheta * c

    # a[0] in Python -> a[1] in Julia
    return a[1] * (be * cf - bphi * sf), a[1] * (be * sf + bphi * cf), a[1] * (br * c - btheta * s)
end

function r_s(a, r, theta)
    # Python: a[1]/r + a[2]*r/sqrt(r^2+a[10]^2) + a[3]*r/(r^2+a[11]^2) + ...
    # Julia indices: Python a[i] -> Julia a[i+1]
    cost, cos2t = cos(theta), cos(2.0 * theta)
    return r + a[2] / r + a[3] * r / sqrt(r^2 + a[11]^2) + a[4] * r / (r^2 + a[12]^2) +
        (a[5] + a[6] / r + a[7] * r / sqrt(r^2 + a[13]^2) + a[8] * r / (r^2 + a[14]^2)) * cost +
        (a[9] * r / sqrt(r^2 + a[15]^2) + a[10] * r / (r^2 + a[16]^2)^2) * cos2t
end

function theta_s(a, r, theta)
    # Python: (a[16] + a[17]/r + a[18]/r^2 + a[19]*r/sqrt(r^2+a[26]^2))*sin(theta) + ...
    # Julia indices: Python a[i] -> Julia a[i+1]
    sint, sin2t, sin3t = sin(theta), sin(2.0 * theta), sin(3.0 * theta)
    return theta + (a[17] + a[18] / r + a[19] / r^2 + a[20] * r / sqrt(r^2 + a[27]^2)) * sint +
        (a[21] + a[22] * r / sqrt(r^2 + a[28]^2) + a[23] * r / (r^2 + a[29]^2)) * sin2t +
        (a[24] + a[25] / r + a[26] * r / (r^2 + a[30]^2)) * sin3t
end

function fialcos(r, theta, phi, n, theta0, dt)
    # Full implementation matching Python fialcos
    sinte = sin(theta)
    ro = r * sinte
    coste = cos(theta)
    sinfi, cosfi = sin(phi), cos(phi)
    tg = sinte / (1.0 + coste)   # tan(theta/2)
    ctg = sinte / (1.0 - coste)  # cot(theta/2)

    tetanp = theta0 + dt
    tetanm = theta0 - dt

    tgp, tgm, tgm2, tgp2 = 0.0, 0.0, 0.0, 0.0
    if theta >= tetanm
        tgp = tan(tetanp * 0.5)
        tgm = tan(tetanm * 0.5)
        tgm2 = tgm * tgm
        tgp2 = tgp * tgp
    end

    cosm1, sinm1 = 1.0, 0.0
    tm = 1.0
    tgm2m, tgp2m = 1.0, 1.0

    btn, bpn = 0.0, 0.0
    ccos_m, ssin_m = 0.0, 0.0

    for m in 1:n
        tm = tm * tg
        ccos_m = cosm1 * cosfi - sinm1 * sinfi
        ssin_m = sinm1 * cosfi + cosm1 * sinfi
        cosm1, sinm1 = ccos_m, ssin_m

        if theta < tetanm
            t = tm
            dtt = 0.5 * m * tm * (tg + ctg)
        elseif theta < tetanp
            tgm2m = tgm2m * tgm2
            fc = 1.0 / (tgp - tgm)
            fc1 = 1.0 / (2.0 * m + 1.0)
            tgm2m1 = tgm2m * tgm
            tg21 = 1.0 + tg * tg
            t = fc * (tm * (tgp - tg) + fc1 * (tm * tg - tgm2m1 / tm))
            dtt = 0.5 * m * fc * tg21 * (tm / tg * (tgp - tg) - fc1 * (tm - tgm2m1 / (tm * tg)))
        else
            tgp2m = tgp2m * tgp2
            tgm2m = tgm2m * tgm2
            fc = 1.0 / (tgp - tgm)
            fc1 = 1.0 / (2.0 * m + 1.0)
            t = fc * fc1 * (tgp2m * tgp - tgm2m * tgm) / tm
            dtt = -t * m * 0.5 * (tg + ctg)
        end

        btn = m * t * ccos_m / ro
        bpn = -dtt * ssin_m / r
    end

    btheta = btn * 800.0
    bphi = bpn * 800.0

    return btheta, bphi
end

function birk_shl(a, ps, x_sc, x, y, z)
    # Python indices: a[84], a[85] -> Julia a[85], a[86]
    cps, sps = cos(ps), sin(ps)
    s3ps = 2.0 * cps
    pst1, pst2 = ps * a[85], ps * a[86]
    st1, ct1 = sin(pst1), cos(pst1)
    st2, ct2 = sin(pst2), cos(pst2)
    x1, z1 = x * ct1 - z * st1, x * st1 + z * ct1
    x2, z2 = x * ct2 - z * st2, x * st2 + z * ct2

    bx, by, bz = 0.0, 0.0, 0.0
    l = 0
    for m in 1:2
        for i in 1:3
            # Python: a[71+i], a[77+i] -> Julia: a[72+i], a[78+i]
            p, q = a[72 + i], a[78 + i]
            cypi, sypi = cos(y / p), sin(y / p)
            cyqi, syqi = cos(y / q), sin(y / q)
            for k in 1:3
                # Python: a[74+k], a[80+k] -> Julia: a[75+k], a[81+k]
                r, s = a[75 + k], a[81 + k]
                szrk, czrk = sin(z1 / r), cos(z1 / r)
                czsk, szsk = cos(z2 / s), sin(z2 / s)
                sqpr = sqrt(1.0 / p^2 + 1.0 / r^2)
                sqqs = sqrt(1.0 / q^2 + 1.0 / s^2)
                epr = exp(x1 * sqpr)
                eqs = exp(x2 * sqqs)
                for n in 1:2
                    for nn in 1:2
                        if m == 1
                            fx = -sqpr * epr * cypi * szrk
                            fy = epr * sypi * szrk / p
                            fz = -epr * cypi * czrk / r
                            if n == 1
                                hx, hy, hz = nn == 1 ? (fx, fy, fz) : (fx * x_sc, fy * x_sc, fz * x_sc)
                            else
                                hx, hy, hz = nn == 1 ? (fx * cps, fy * cps, fz * cps) : (fx * cps * x_sc, fy * cps * x_sc, fz * cps * x_sc)
                            end
                        else  # m == 2
                            fx = -sps * sqqs * eqs * cyqi * czsk
                            fy = sps / q * eqs * syqi * czsk
                            fz = sps / s * eqs * cyqi * szsk
                            if n == 1
                                hx, hy, hz = nn == 1 ? (fx, fy, fz) : (fx * x_sc, fy * x_sc, fz * x_sc)
                            else
                                hx, hy, hz = nn == 1 ? (fx * s3ps, fy * s3ps, fz * s3ps) : (fx * s3ps * x_sc, fy * s3ps * x_sc, fz * s3ps * x_sc)
                            end
                        end
                        l += 1
                        ct, st = m == 1 ? (ct1, st1) : (ct2, st2)
                        hxr = hx * ct + hz * st
                        hzr = -hx * st + hz * ct
                        bx += hxr * a[l]
                        by += hy * a[l]
                        bz += hzr * a[l]
                    end
                end
            end
        end
    end
    return bx, by, bz
end
