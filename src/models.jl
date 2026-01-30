"""
    (m::TsyganenkoModel)(x, y, z, ps, args...; kw...) -> (Bx, By, Bz)

Compute GSM components of the external magnetic field [nT] using Tsyganenko model `m`, given position in GSM coordinates (`x, y, z`) [Earth radii] and geodipole tilt angle [radians] `ps`.
"""
abstract type TsyganenkoModel <: ExternalFieldModel end

getcsys(::TsyganenkoModel) = (GSM(), Cartesian3())

# Generic model types
"""
    T89(iopt::Integer)

Tsyganenko 89 model with Kp-dependent parameters.

# Parameters
- `iopt`: Ground disturbance level index (1-7)
  - 1: Kp = 0, 0+
  - 2: Kp = 1-, 1, 1+
  - 3: Kp = 2-, 2, 2+
  - 4: Kp = 3-, 3, 3+
  - 5: Kp = 4-, 4, 4+
  - 6: Kp = 5-, 5, 5+
  - 7: Kp ≥ 6-

# Model Description
Valid up to geocentric distances of 70 RE. Based on merged IMP-A through J (1966-1974),
HEOS-1 and -2 (1969-1974), and ISEE-1 and -2 spacecraft data.

# Usage
```julia
model = T89(2)  # Kp level 2
B = model([1, 2, 3], ps)  # Compute field at position
```

# References
- Tsyganenko, N.A., "A magnetospheric magnetic field model with a warped tail current sheet",
  Planet. Space Sci., 37, 5-20, 1989.
- [Fortran implementation](https://geo.phys.spbu.ru/~tsyganenko/models/t89/T89d_dp.for)
"""
struct T89{T <: Integer, C} <: TsyganenkoModel
    iopt::T
    cache::C
end

function T89(iopt)
    @assert 1 ≤ iopt ≤ 7 "iopt must be between 1 and 7"
    cache = t89_init(iopt)
    return T89(iopt, cache)
end

T89(; iopt) = T89(iopt)

"""
    T96(; pdyn, dst, byimf, bzimf)

Tsyganenko 96 model with solar wind parameters.

# Parameters
- `pdyn`: Solar wind dynamic pressure [nPa]
- `dst`: Dst index [nT]
- `byimf`: IMF By component [nT]
- `bzimf`: IMF Bz component [nT]

# Usage
```julia
model = T96(pdyn=2.0, dst=-87.0, byimf=2.0, bzimf=-5.0)
B = model([1, 2, 3], ps)  # Compute field at position
```

# Model Description
Data-based model calibrated by solar wind pressure, Dst index, and IMF By/Bz components.
Includes realistic magnetopause, Region 1 and 2 Birkeland current systems, and IMF penetration.

Valid parameter ranges (caution needed outside these ranges):
- Pdyn: 0.5 to 10 nPa
- Dst: -100 to +20 nT
- ByIMF and BzIMF: -10 to +10 nT

# References
- Tsyganenko & Stern, JGR, v.101, p.27187-27198, 1996
- Tsyganenko, JGR, v.100, p.5599-5612, 1995
"""
struct T96{T} <: TsyganenkoModel
    pdyn::T   # Solar wind dynamic pressure [nPa]
    dst::T    # Dst index [nT]
    byimf::T  # IMF By [nT]
    bzimf::T  # IMF Bz [nT]
end

function T96(; pdyn, dst, byimf, bzimf)
    return T96(promote(pdyn, dst, byimf, bzimf)...)
end

"""
    T01(; pdyn, dst, byimf, bzimf, g1=0.0, g2=0.0)

Tsyganenko 01 model with solar wind parameters.

# Parameters
- `pdyn`: Solar wind dynamic pressure [nPa]
- `dst`: Dst index [nT]
- `byimf`: IMF By component [nT]
- `bzimf`: IMF Bz component [nT]
- `g1`: G1 index (default: 0.0)
- `g2`: G2 index (default: 0.0)

# Usage
```julia
model = T01(pdyn=2.0, dst=-87.0, byimf=2.0, bzimf=-5.0)
B = model([1, 2, 3], ps)  # Compute field at position
```

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
struct T01{T} <: TsyganenkoModel
    pdyn::T
    dst::T
    byimf::T
    bzimf::T
    g1::T
    g2::T
end

function T01(; pdyn, dst, byimf, bzimf, g1 = 0.0, g2 = 0.0)
    return T01(promote(pdyn, dst, byimf, bzimf, g1, g2)...)
end
