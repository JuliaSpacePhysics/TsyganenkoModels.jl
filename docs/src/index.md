```@meta
CurrentModule = TsyganenkoModels
```

# TsyganenkoModels

Modeling of Earth's Magnetosphere Using Spacecraft Magnetometer Data.

## Features

- Magnetic field model (native Julia implementations)
  - [x] T89: A magnetospheric magnetic field model with a warped tail current sheet
  - [x] T96: Effects of the solar wind conditions on the global magnetospheric configuration
  - [x] T01/T02: A model of the near magnetosphere with a dawn-dusk asymmetry
  - [x] TS05/TS04: a dynamical empirical model of the inner storm-time magnetosphere

!!! note "IRBEM.jl"
    [`IRBEM.jl`](https://github.com/JuliaSpacePhysics/IRBEM.jl) is a Julia wrapper for the IRBEM Fortran library that exposes magnetic field computation via `GET_FIELD_MULTI` and supports more Tsyganenko models, but may be outdated and slower than the native Julia implementations.

## Installation

```julia
using Pkg
Pkg.add("TsyganenkoModels")
```

## Quickstart

```@example quickstart
using TsyganenkoModels
using Dates

t = DateTime("1970-01-01T00:01:40")
𝐫 = [1, 2, 3]
iopt = 2
db_jl_t = t89(𝐫, t, iopt)

pdyn = 2.0   # Solar wind dynamic pressure [nPa]
dst = -87.0  # Dst index [nT]
byimf = 2.0  # IMF By [nT]
bzimf = -5.0 # IMF Bz [nT]
ps = -0.533585131
result_t01 = t01(𝐫, ps, pdyn, dst, byimf, bzimf)
result_t01 = t96(𝐫, ps, pdyn, dst, byimf, bzimf)
```

## Comparison with [`geopack`](https://github.com/tsssss/geopack)

```@example comparison
using TsyganenkoModels
using GeoCotrans
using Test
using Geopack
using Dates

t = DateTime("1970-01-01T00:01:40")
𝐫 = [1, 2, 3]

# using geopack
xgsm, ygsm, zgsm = 𝐫
ut = datetime2unix(t)    # 1970-01-01/00:01:40 UT.
ps = Geopack.recalc(ut)
db_py = Geopack.t89(2, ps, xgsm, ygsm, zgsm)

# using TsyganenkoModels
db_jl = t89(𝐫, ps, 2)

@test db_py == db_jl
```

### Model Result Comparison

```@example comparison
# T01 model comparison
pdyn = 2.0   # Solar wind dynamic pressure [nPa]
dst = -87.0  # Dst index [nT]
byimf = 2.0  # IMF By [nT]
bzimf = -5.0 # IMF Bz [nT]
g1 = 0.0
g2 = 0.0
parmod = [pdyn, dst, byimf, bzimf, g1, g2]

db_t96_py = Geopack.t96(parmod, ps, xgsm, ygsm, zgsm)
db_t96_jl = t96(𝐫, ps, pdyn, dst, byimf, bzimf)
@test collect(db_t96_jl) ≈ db_t96_py rtol = 1e-6

db_t01_py = Geopack.t01(parmod, ps, xgsm, ygsm, zgsm)
db_t01_jl = t01(𝐫, ps, pdyn, dst, byimf, bzimf)

@test collect(db_t01_jl) ≈ db_t01_py rtol = 1e-6
```

The internal magnetic field could be computed using `igrf` from [`GeoCotrans`](https://github.com/JuliaSpacePhysics/GeoCotrans.jl).

Here we verify the result with `igrf_gsm` from `Geopack`.

```@repl comparison
using GeoCotrans

b0_py = Geopack.igrf_gsm(xgsm, ygsm, zgsm)
b0_jl = GeoCotrans.igrf(GSM(𝐫) .* GeoCotrans.R🜨, t)
@test b0_jl ≈ b0_py rtol = 1e-5
@test GeoCotrans.get_igrf_coeffs(t)[1] ≈ Geopack.load_igrf(ut)[1]
```

### Performance Benchmarks

```@repl comparison
using Chairmarks
@b TsyganenkoModels.t89(𝐫, ps, 2), Geopack.t89(2, ps, xgsm, ygsm, zgsm)
@b TsyganenkoModels.t96(𝐫, ps, pdyn, dst, byimf, bzimf), Geopack.t96([pdyn, dst, byimf, bzimf, 0, 0], ps, xgsm, ygsm, zgsm)
@b TsyganenkoModels.t01(𝐫, ps, pdyn, dst, byimf, bzimf), Geopack.t01([pdyn, dst, byimf, bzimf, 0, 0], ps, xgsm, ygsm, zgsm)
@b GeoCotrans.igrf(GSM(xgsm, ygsm, zgsm) .* GeoCotrans.R🜨, t), Geopack.igrf_gsm(xgsm, ygsm, zgsm)
```

## API

```@index
```

```@autodocs
Modules = [TsyganenkoModels]
```