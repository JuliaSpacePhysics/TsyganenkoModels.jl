# TsyganenkoModels

[![Build Status](https://github.com/JuliaSpacePhysics/TsyganenkoModels.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JuliaSpacePhysics/TsyganenkoModels.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/JuliaSpacePhysics/TsyganenkoModels.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaSpacePhysics/TsyganenkoModels.jl)
[![Aqua](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

Modeling of Earth's Magnetosphere Using Spacecraft Magnetometer Data

**Installation**: at the Julia REPL, run `using Pkg; Pkg.add("TsyganenkoModels")`

**Documentation**: [![Dev](https://img.shields.io/badge/docs-dev-blue.svg?logo=julia)](https://JuliaSpacePhysics.github.io/TsyganenkoModels.jl/dev/)

## Features and Roadmap

- Magnetic field model
  - [x] Supported models: T89, T96 (native Julia implementations).
  - [x] T02: A model of the near magnetosphere with a dawn-dusk asymmetry
  - [ ] TA15 model: A forecasting model of the magnetosphere driven by optimal solar-wind coupling functions
  - [ ] TA16 model: An empirical RBF model of the magnetosphere parameterized by interplanetary and ground-based drivers
- Plasma model
  - [ ] Tsyganenko and Mukai (2003): a simple analytical model of the central plasma sheet ion parameters (10 − 50 RE), based on Geotail data.

## Elsewhere

- [tsyganenko/empirical-models](https://geo.phys.spbu.ru/~tsyganenko/empirical-models/): original source of the Tsyganenko models in Fortran
- [tsssss/geopack](https://github.com/tsssss/geopack): Python version of geopack and Tsyganenko models

A Julia wrapper for [`geopack`](https://github.com/tsssss/geopack) is available in the `lib/Geopack.jl` directory and can be installed with `using Pkg; Pkg.develop(url="https://github.com/JuliaSpacePhysics/TsyganenkoModels.jl", subdir="lib/Geopack.jl")`. It is mainly used for testing and benchmarking.