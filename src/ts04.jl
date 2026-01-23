# Tsyganenko-Sitnov 04 Magnetospheric Magnetic Field Model (TS04)

module TS04Impl
    using ..TsyganenkoModels: birk_tot, full_rc, dipole, deformed, shlcar3x3
    include("ts04_consts.jl")
    include("ts04_funcs.jl")
end

using .TS04Impl

"""
    ts04(x, y, z, ps, (pdyn, dst, byimf, bzimf, w1, w2, w3, w4, w5, w6)) -> (Bx, By, Bz)

Compute GSM components of the external magnetic field [nT] using the Tsyganenko-Sitnov 04 model, given position in GSM coordinates (`x, y, z`) [Earth radii], geodipole tilt angle [radians] `ps`, and model parameters.
- `pdyn`: Solar wind dynamic pressure [nPa]
- `dst`: Dst index [nT]
- `byimf`: IMF By component [nT]
- `bzimf`: IMF Bz component [nT]
- `w1-w6`: Time integrals of the driving variables (time integrals from the beginning of a storm, see reference for definitions)

# References
- Tsyganenko, N. A., and Sitnov, M. I. (2005), Modeling the dynamics of the inner magnetosphere during strong geomagnetic storms, J. Geophys. Res., 110, A03208, https://doi.org/10.1029/2004JA010798.
- [TS05 magnetic field model](https://geo.phys.spbu.ru/~tsyganenko/empirical-models/magnetic_field/ts05)
- [Fortran implementation (TS04c.for)](https://geo.phys.spbu.ru/~tsyganenko/models/ts04/TS04c.for)
"""
ts04(x, y, z, ps, params) = _ts04(x, y, z, ps, params...)

function _ts04(x, y, z, ps, pdyn, dst, byimf, bzimf, w1 = 0.0, w2 = 0.0, w3 = 0.0, w4 = 0.0, w5 = 0.0, w6 = 0.0)
    if x < -20.0
        @warn "The model is valid sunward from X=-15 Re only, while you are trying to use it at X=$x"
    end
    dst_ast = dst * 0.8 - 13.0 * sqrt(pdyn)
    return TS04Impl.extall(pdyn, dst_ast, byimf, bzimf, w1, w2, w3, w4, w5, w6, ps, x, y, z)
end
