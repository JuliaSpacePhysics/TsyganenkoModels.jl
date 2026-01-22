```@meta
CurrentModule = TsyganenkoModels
```

# TsyganenkoModels

Documentation for [TsyganenkoModels](https://github.com/JuliaSpacePhysics/TsyganenkoModels.jl).

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

### T96/T01 Model Comparison

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
db_t01_jl = t01(𝐫, ps, pdyn, dst, byimf, bzimf; g1, g2)

@test collect(db_t01_jl) ≈ db_t01_py rtol = 1e-6
```

The internal magnetic field could be computed using `igrf_Bgsm` from [`GeoCotrans`](https://github.com/JuliaSpacePhysics/GeoCotrans.jl).

Here we verify the result with `igrf_gsm` from `Geopack`.

```@repl comparison
using GeoCotrans

b0_py = Geopack.igrf_gsm(xgsm, ygsm, zgsm)
b0_jl = GeoCotrans.igrf_B(GSM(𝐫) .* GeoCotrans.R🜨, t)
@test b0_jl ≈ b0_py rtol = 1e-5
@test GeoCotrans.get_igrf_coeffs(t)[1] ≈ Geopack.load_igrf(ut)[1]
```

Benchmarks

```@repl comparison
using Chairmarks
@b TsyganenkoModels.t89(𝐫, ps, 2), Geopack.t89(2, ps, xgsm, ygsm, zgsm)
@b GeoCotrans.igrf_B(GSM(xgsm, ygsm, zgsm) .* GeoCotrans.R🜨, t), Geopack.igrf_gsm(xgsm, ygsm, zgsm)
```

## API

```@index
```

```@autodocs
Modules = [TsyganenkoModels]
```