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