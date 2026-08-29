# TsyganenkoModels

Modeling of Earth's Magnetosphere Using Spacecraft Magnetometer Data

## Quickstart

```julia
using Pkg; Pkg.add("TsyganenkoModels")
using TsyganenkoModels
using Dates

# Create model configurations
model_t89 = T89(2)  # Kp level 2
param = (; pdyn=2.0, dst=-87.0, byimf=2.0, bzimf=-5.0)
# pdyn: Solar wind dynamic pressure [nPa]

# Calculate fields at position
t = DateTime("1970-01-01T00:01:40")
𝐫 = [1, 2, 3]
ps = -0.533585131  # dipole tilt angle [radians]

B_t89 = model_t89(𝐫, ps)
B_t96 = T96(param)(𝐫, ps)
B_t01 = T01(param)(𝐫, ps)
B_ts04 = TS04(param)(𝐫, ps)

# Using `TsyIGRF` to combine IGRF14 model and Tsyganenko model, default TsyIGRF() == TsyIGRF(T89(iopt=3))
TsyIGRF()(𝐫, t)
```

## Features and Roadmap

- Magnetic field model
  - [x] T89: A magnetospheric magnetic field model with a warped tail current sheet
  - [x] T96: Effects of the solar wind conditions on the global magnetospheric configuration
  - [x] T01/T02: A model of the near magnetosphere with a dawn-dusk asymmetry
  - [x] TS05/TS04: a dynamical empirical model of the inner storm-time magnetosphere
  - [ ] TA15 model: A forecasting model of the magnetosphere driven by optimal solar-wind coupling functions
  - [ ] TA16 model: An empirical RBF model of the magnetosphere parameterized by interplanetary and ground-based drivers
- Plasma model
  - [ ] Tsyganenko and Mukai (2003): a simple analytical model of the central plasma sheet ion parameters (10 − 50 RE), based on Geotail data.

For Python users, see the wrapper in [`python/`](python/README.md) (PyPI: [`tsyganenkomodels-jl`](https://pypi.org/project/tsyganenkomodels-jl/)).

## Elsewhere

- [tsyganenko/empirical-models](https://geo.phys.spbu.ru/~tsyganenko/empirical-models/): original source of the Tsyganenko models in Fortran
- [tsssss/geopack](https://github.com/tsssss/geopack): Python version of geopack and Tsyganenko models
- [IDL Geopack DLM](https://korthhaus.com/idl-software/idl-geopack-dlm/)

A Julia wrapper for [`geopack`](https://github.com/tsssss/geopack) is available in the `lib/Geopack.jl` directory and can be installed with `using Pkg; Pkg.develop(url="https://github.com/JuliaSpacePhysics/TsyganenkoModels.jl", subdir="lib/Geopack.jl")`. It is mainly used for testing and benchmarking.

## Status

[![DOI](https://zenodo.org/badge/1139356853.svg)](https://doi.org/10.5281/zenodo.18435633)
[![version](https://juliahub.com/docs/General/TsyganenkoModels/stable/version.svg)](https://juliahub.com/ui/Packages/General/TsyganenkoModels)
![PyPI - Version](https://img.shields.io/pypi/v/tsyganenkomodels-jl)
[![Coverage](https://codecov.io/gh/JuliaSpacePhysics/TsyganenkoModels.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaSpacePhysics/TsyganenkoModels.jl)