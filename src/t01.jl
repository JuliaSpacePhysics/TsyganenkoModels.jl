#
# Tsyganenko 01 Magnetospheric Magnetic Field Model (T01_01)
#
# Based on the Fortran implementation by N.A. Tsyganenko
# https://geo.phys.spbu.ru/~tsyganenko/models/t01/T01_01c.for
#
# References:
# - Tsyganenko, N.A., "A new data-based model of the near magnetosphere magnetic field:
#   1. Mathematical structure. 2. Parameterization and fitting to observations."
#   J.Geophys.Res., 2002.
#

module T01Impl
    using ..TsyganenkoModels: birk_tot, full_rc, dipole, deformed, shlcar3x3
    include("t01_consts.jl")
    include("t01_funcs.jl")
end

using .T01Impl

"""
    t01(x, y, z, ps, pdyn, dst, byimf, bzimf; g1 = 0.0, g2 = 0.0) -> (Bx, By, Bz)

Compute GSM components of the external magnetic field using the Tsyganenko 01 model.

# Parameters
- `x, y, z`: Position in GSM coordinates [Earth radii]
- `ps`: Geodipole tilt angle [radians]
- `pdyn`: Solar wind dynamic pressure [nPa]
- `dst`: Dst index [nT]
- `byimf`: IMF By component [nT]
- `bzimf`: IMF Bz component [nT]
- `g1`: G1 index (see Tsyganenko 2001 for definition)
- `g2`: G2 index (see Tsyganenko 2001 for definition)

# Returns
- `(Bx, By, Bz)`: Magnetic field components in GSM coordinates [nT]

# Model Description
Data-based model of the external magnetospheric magnetic field, calibrated by:
1. Solar wind pressure PDYN (nanopascals)
2. Dst (nanotesla)
3. IMF By and Bz (nanotesla)
4. G1 and G2 indices

**ATTENTION**: The model is based on data taken sunward from X=-15 Re, and hence
becomes invalid at larger tailward distances.

# References
- Tsyganenko, N.A., JGR, 2002
"""
function t01(x, y, z, ps, pdyn, dst, byimf, bzimf; g1 = 0.0, g2 = 0.0)
    if x < -20.0
        @warn "The model is valid sunward from X=-15 Re only, while you are trying to use it at X=$x"
    end
    dst_ast = dst * 0.8 - 13.0 * sqrt(pdyn)
    return T01Impl.extall(pdyn, dst_ast, byimf, bzimf, g1, g2, ps, x, y, z)
end

function t01(x, y, z, t::AbstractTime, args...; kw...)
    ps = dipole_tilt(t)
    return t01(x, y, z, ps, args...; kw...)
end

@inline function t01(r, args...; kw...)
    @assert length(r) == 3
    return t01(r[1], r[2], r[3], args...; kw...)
end
