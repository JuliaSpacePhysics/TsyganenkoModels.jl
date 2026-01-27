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
    t01(x, y, z, ps, pdyn, dst, byimf, bzimf, g1 = 0.0, g2 = 0.0) -> (Bx, By, Bz)

Compute GSM components of the external magnetic field using the Tsyganenko 01 model.
"""
function t01(x, y, z, ps, pdyn, dst, byimf, bzimf, g1 = 0.0, g2 = 0.0)
    x < -20.0 && @warn "The model is valid sunward from X=-15 Re only, while you are trying to use it at X=$x"
    dst_ast = dst * 0.8 - 13.0 * sqrt(pdyn)
    return GSM(T01Impl.extall(pdyn, dst_ast, byimf, bzimf, g1, g2, ps, x, y, z))
end

# Functor interfaces for T01
evalmodel(m::T01, x, y, z, ps) = t01(x, y, z, ps, m.pdyn, m.dst, m.byimf, m.bzimf, m.g1, m.g2)
