# TsyganenkoModels

Python wrapper for [TsyganenkoModels.jl](https://github.com/JuliaSpacePhysics/TsyganenkoModels.jl) - evaluate Tsyganenko external magnetic field models from Python.

## Installation

```bash
pip install tsyganenkomodels-jl
```

The Julia package installs automatically on first use.

## Usage

```python
import tsyganenkomodels as tm
import datetime

# Compute the dipole tilt angle (radians)
ps = tm.dipole_tilt(datetime.datetime(2001, 1, 1, 2, 3, 4))

# Instantiate a model and evaluate GSM field (nT)
model = tm.T96(pdyn=2.0, dst=-87.0, byimf=2.0, bzimf=-5.0)
b = model([-5.1, 0.3, 2.8], ps)
```

## API

Available constructors/functions mirror the Julia package:

- `T89(iopt)`
- `T96(pdyn, dst, byimf, bzimf)`
- `T01(pdyn, dst, byimf, bzimf, g1=0.0, g2=0.0)`
- `TS04(pdyn, dst, byimf, bzimf, g1, g2, g3, w1, w2, w3, w4, w5, w6)`
- `t89`, `t96`, `t01`, `ts04` (functional wrappers)
- `dipole_tilt(t)`

See the [TsyganenkoModels.jl documentation](https://JuliaSpacePhysics.github.io/TsyganenkoModels.jl) for full details.
