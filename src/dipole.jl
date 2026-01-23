"""
    dipole_tilt(time) -> ψ

Compute the dipole tilt angle ψ (radians): ψ = asin(ẑ_dipole · x̂_sun).
"""
function dipole_tilt(t)
    _, ra, dec, _, _ = GeoCotrans.csundir(t)
    sun_gei = GeoCotrans.calc_sun_gei(ra, dec)
    dipole_gei = GeoCotrans.calc_dipole_gei(t)
    sps = clamp(dot(dipole_gei, sun_gei), -1.0, 1.0)
    return asin(sps)
end

function dipdistr(x, y, z, mode)
    x2 = x^2
    rho2 = x2 + y^2
    r2 = rho2 + z^2
    r3 = r2 * sqrt(r2)

    if mode == 0
        bx = z / rho2^2 * (r2 * (y^2 - x2) - rho2 * x2) / r3
        by = -x * y * z / rho2^2 * (2 * r2 + rho2) / r3
        bz = x / r3
    else
        bx = z / rho2^2 * (y^2 - x2)
        by = -2 * x * y * z / rho2^2
        bz = x / rho2
    end
    return bx, by, bz
end

function dipole(ps, x, y, z; q0 = 30574)
    sps, cps = sincos(ps)

    p = x^2
    u = z^2
    v = 3 * z * x
    t = y^2
    q = q0 / sqrt(p + t + u)^5
    bx = q * ((t + u - 2 * p) * sps - v * cps)
    by = -3 * y * q * (x * sps + z * cps)
    bz = q * ((p + t - 2 * u) * cps - v * sps)

    return bx, by, bz
end

# This subroutine returns the shielding field for the earth's dipole, represented by
# 2x3x3=18 "cartesian" harmonics, tilted with respect to the z=0 plane
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

# shlcar3x3 coefficients (50 elements)
const SHLCAR_A = [
    -901.2327248, 895.8011176, 817.6208321, -845.5880889, -83.73539535,
    86.58542841, 336.8781402, -329.3619944, -311.294712, 308.6011161,
    31.94469304, -31.30824526, 125.8739681, -372.3384278, -235.4720434,
    286.7594095, 21.86305585, -27.42344605, -150.4874688, 2.669338538,
    1.395023949, -0.5540427503, -56.85224007, 3.681827033, -43.48705106,
    5.103131905, 1.073551279, -0.6673083508, 12.21404266, 4.177465543,
    5.799964188, -0.3977802319, -1.044652977, 0.570356001, 3.536082962,
    -3.222069852, 9.620648151, 6.082014949, 27.75216226, 12.44199571,
    5.122226936, 6.982039615, 20.12149582, 6.150973118, 4.663639687,
    15.73319647, 2.303504968, 5.840511214, 0.08385953499, 0.3477844929,
]