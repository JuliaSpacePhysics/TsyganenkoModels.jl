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

The model interface allows you to configure a model once and use it for multiple field calculations:

```@example quickstart
using TsyganenkoModels
using Dates

# Create model configurations
model_t89 = T89(2)  # Kp level 2
param = (; pdyn=2.0, dst=-87.0, byimf=2.0, bzimf=-5.0)
# pdyn: Solar wind dynamic pressure [nPa]
# dst: Dst index [nT]
# byimf: IMF By [nT]
# bzimf: IMF Bz [nT]

# Calculate fields at position
t = DateTime("1970-01-01T00:01:40")
𝐫 = [1, 2, 3]
ps = -0.533585131  # dipole tilt angle [radians]

# Using dipole tilt angle
B_t89 = T89(2)(𝐫, ps)
```

```@repl quickstart
# Compare with other models
B_t96 = T96(param)(𝐫, ps)
B_t01 = T01(param)(𝐫, ps)
B_ts04 = TS04(param)(𝐫, ps)

# Using time (auto-calculates dipole tilt)
T89(2)(𝐫, t)
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
opt = 2
ut = datetime2unix(t)
ps = Geopack.recalc(ut)
db_py = Geopack.t89(opt, ps, xgsm, ygsm, zgsm)

# using TsyganenkoModels
db_jl = T89(opt)(𝐫, ps)

@test db_py ≈ db_jl
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
db_t96_jl = T96(; pdyn, dst, byimf, bzimf)(𝐫, ps)
@test db_t96_jl ≈ db_t96_py rtol = 1e-6

db_t01_py = Geopack.t01(parmod, ps, xgsm, ygsm, zgsm)
db_t01_jl = T01(; pdyn, dst, byimf, bzimf)(𝐫, ps)

@test db_t01_jl ≈ db_t01_py rtol = 1e-6
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
@b TsyganenkoModels.T89(2)($𝐫, $ps), Geopack.t89(2, $ps, $𝐫...)
@b TsyganenkoModels.T96($pdyn, $dst, $byimf, $bzimf)($𝐫, $ps), Geopack.t96([$pdyn, $dst, $byimf, $bzimf, 0, 0], $ps, $𝐫...)
@b TsyganenkoModels.T01(; pdyn, dst, byimf, bzimf)($𝐫, $ps), Geopack.t01([pdyn, dst, byimf, bzimf, 0, 0], $ps, $𝐫...)

@b GeoCotrans.igrf(GSM($𝐫) .* GeoCotrans.R🜨, t), Geopack.igrf_gsm($𝐫...)
```

## API

```@index
```

```@autodocs
Modules = [TsyganenkoModels]
```