# Reference equatorial dipole field at 1 Earth radius, in nT, used by T96.
const T96_DIPOLE_STRENGTH = 30574

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

function dipole(ps, x, y, z; q0 = T96_DIPOLE_STRENGTH)
    sps, cps = sincos(ps)
    return dipxyz(x, y, z, (-sps, 0, -cps); q0)
end

@inline function dipxyz(x, y, z, m; q0 = T96_DIPOLE_STRENGTH)
    mx, my, mz = m
    r2 = x^2 + y^2 + z^2
    scale = q0 / (r2^2 * sqrt(r2))
    rm = 3 * (x * mx + y * my + z * mz)
    return scale * (x * rm - r2 * mx),
        scale * (y * rm - r2 * my),
        scale * (z * rm - r2 * mz)
end


# This subroutine returns the shielding field for the earth's dipole, represented by
# 2x3x3=18 "cartesian" harmonics, tilted with respect to the z=0 plane
function dipole_shield(x, y, z, ps)
    a = SHLCAR_A
    sps, cps = sincos(ps)
    p, r, q, s = (ntuple(i -> a[36 + 3j + i], 3) for j in 0:3)
    c1(i, k) = a[6i + 2k - 7] + a[6i + 2k - 6] * cps
    c2(i, k) = a[6i + 2k + 11] + a[6i + 2k + 12] * 2cps
    b1 = tilted(ps * a[49], x, y, z) do x, y, z
        cartharm(c1, p, r[1:2], x, y, z, false) .+ cartharm_dr(i -> c1(i, 3), p, r[3], x, y, z)
    end
    b2 = tilted(ps * a[50], x, y, z) do x, y, z
        cartharm(c2, q, s, x, y, z, true)
    end
    return b1 .+ sps .* b2
end

# Perpendicular harmonics of `cartharm` replaced by -r²∂U/∂r (still harmonic: ∇²U = 0 for every
# r); the dipole shield uses them for its third z-wavelength.
function cartharm_dr(c, p, r, x, y, z)
    sz, cz = sincos(z / r)
    bx = by = bz = 0.0
    for i in eachindex(p)
        sy, cy = sincos(y / p[i])
        sq = sqrt(1 / p[i]^2 + 1 / r^2)
        e = c(i) * exp(x * sq)
        bx -= e * cy * (sq * z * cz + sz / r * (x + 1 / sq))
        by += e / p[i] * sy * (z * cz + x / r * sz / sq)
        bz -= e * cy * (cz * (1 + x / r^2 / sq) - z / r * sz)
    end
    return bx, by, bz
end

# dipole_shield coefficients (50 elements)
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
