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

function dipole(ps, x, y, z)
    sps, cps = sincos(ps)

    p = x^2
    u = z^2
    v = 3 * z * x
    t = y^2
    q = 30574 / sqrt(p + t + u)^5
    bx = q * ((t + u - 2 * p) * sps - v * cps)
    by = -3 * y * q * (x * sps + z * cps)
    bz = q * ((p + t - 2 * u) * cps - v * sps)

    return bx, by, bz
end