# Tsyganenko-Sitnov 04 Magnetospheric Magnetic Field Model (TS04)

"""
    TS04(; pdyn, dst, byimf, bzimf, w1=0.0, w2=0.0, w3=0.0, w4=0.0, w5=0.0, w6=0.0)

Tsyganenko-Sitnov 04 model with solar wind parameters and W indices.

# Parameters
- `pdyn`: Solar wind dynamic pressure [nPa]
- `dst`: Dst index [nT]
- `byimf`: IMF By component [nT]
- `bzimf`: IMF Bz component [nT]
- `w1-w6`: Time integrals of the driving variables (time integrals from the beginning of a storm, see reference for definitions)

# Usage
```julia
model = TS04(pdyn=2.0, dst=-87.0, byimf=2.0, bzimf=-5.0, w1=0.5)
B = model([1, 2, 3], ps)  # Compute field at position
```

# References
- Tsyganenko, N. A., and Sitnov, M. I. (2005), Modeling the dynamics of the inner magnetosphere during strong geomagnetic storms, J. Geophys. Res., 110, A03208, https://doi.org/10.1029/2004JA010798.
- [TS05 magnetic field model](https://geo.phys.spbu.ru/~tsyganenko/empirical-models/magnetic_field/ts05)
- [Fortran implementation (TS04c.for)](https://geo.phys.spbu.ru/~tsyganenko/models/ts04/TS04c.for)
"""
struct TS04{T} <: TsyganenkoModel
    pdyn::T
    dst::T
    byimf::T
    bzimf::T
    w1::T
    w2::T
    w3::T
    w4::T
    w5::T
    w6::T
end

function TS04(; pdyn, dst, byimf, bzimf, w1 = 0.0, w2 = 0.0, w3 = 0.0, w4 = 0.0, w5 = 0.0, w6 = 0.0)
    return TS04(promote(pdyn, dst, byimf, bzimf, w1, w2, w3, w4, w5, w6)...)
end

module TS04Impl
    using ..TsyganenkoModels: birk_tot, full_rc, dipole, deformed, shlcar3x3
    include("ts04_consts.jl")
    include("ts04_funcs.jl")
end

using .TS04Impl

function ts04(x, y, z, ps, pdyn, dst, byimf, bzimf, w1 = 0.0, w2 = 0.0, w3 = 0.0, w4 = 0.0, w5 = 0.0, w6 = 0.0)
    if x < -20.0
        @warn "The model is valid sunward from X=-15 Re only, while you are trying to use it at X=$x"
    end
    dst_ast = dst * 0.8 - 13.0 * sqrt(pdyn)
    return GSM(TS04Impl.extall(pdyn, dst_ast, byimf, bzimf, w1, w2, w3, w4, w5, w6, ps, x, y, z))
end

(m::TS04)(x, y, z, ps) = ts04(x, y, z, ps, m.pdyn, m.dst, m.byimf, m.bzimf, m.w1, m.w2, m.w3, m.w4, m.w5, m.w6)
