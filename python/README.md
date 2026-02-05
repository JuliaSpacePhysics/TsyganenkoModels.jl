# TsyganenkoModels Python Wrapper

A Python wrapper for [TsyganenkoModels.jl](https://github.com/JuliaSpacePhysics/TsyganenkoModels.jl), providing implementations of Tsyganenko empirical magnetospheric magnetic field models.

## Installation

```bash
pip install tsyganenko
```

Or install from source:

```bash
cd python
pip install -e .
```

### Prerequisites

This package requires Julia to be installed. The first import will automatically install the required Julia packages.

## Supported Models

- **T89**: Tsyganenko 1989 model (Kp-dependent)
- **T96**: Tsyganenko 1996 model (solar wind dependent)
- **T01**: Tsyganenko 2001 model (with dawn-dusk asymmetry)
- **TS04**: Tsyganenko-Sitnov 2004 model (storm-time dynamics)

## Usage

### Functional API

```python
import numpy as np
from tsyganenko import t89, t96, t01, ts04, dipole_tilt
from datetime import datetime

# Calculate dipole tilt angle
ps = dipole_tilt(datetime(2001, 1, 1, 2, 3, 4))

# Single position calculation
position = [-5.1, 0.3, 2.8]  # GSM coordinates in Earth radii
B = t89(position, ps=ps, iopt=2)
print(f"B = {B}")  # [Bx, By, Bz] in nT

# Array of positions (vectorized)
positions = np.array([
    [-5.1, 0.3, 2.8],
    [-4.0, 0.5, 3.0],
    [-6.0, 0.0, 2.0],
])
B = t89(positions, ps=ps, iopt=2)
print(f"B shape: {B.shape}")  # (3, 3)

# T96 model with solar wind parameters
B = t96(position, ps=ps, pdyn=2.0, dst=-87.0, byimf=2.0, bzimf=-5.0)

# T01 model with G-indices
B = t01(position, ps=ps, pdyn=2.0, dst=-87.0, byimf=2.0, bzimf=-5.0, g1=1.0, g2=0.5)

# TS04 model with W-indices
B = ts04(position, ps=ps, pdyn=2.0, dst=-87.0, byimf=2.0, bzimf=-5.0, w1=0.5, w2=0.3)
```

### Object-Oriented API

```python
from tsyganenko import T89, T96, T01, TS04

# Create model instances
model = T89(iopt=2)
B = model([-5.1, 0.3, 2.8], ps=-0.533)

# Models with solar wind parameters
model = T96(pdyn=2.0, dst=-87.0, byimf=2.0, bzimf=-5.0)
B = model([-5.1, 0.3, 2.8], ps=-0.533)

# Array input also works
positions = np.array([[-5.1, 0.3, 2.8], [-4.0, 0.5, 3.0]])
B = model(positions, ps=-0.533)  # Shape: (2, 3)
```

## Input Formats

The wrapper accepts various Python array-like inputs:

- Python lists: `[x, y, z]` or `[[x1, y1, z1], [x2, y2, z2], ...]`
- Python tuples: `(x, y, z)`
- NumPy arrays: `np.array([x, y, z])` or `np.array([[x1, y1, z1], ...])`

All inputs are converted to float64 for computation.

## Coordinate System

All positions and magnetic field components are in **GSM (Geocentric Solar Magnetospheric)** coordinates:
- Positions in Earth radii (RE)
- Magnetic field in nanotesla (nT)

## Model Parameters

### T89
- `iopt`: Kp index level (1-7)
  - 1: Kp = 0, 0+
  - 2: Kp = 1-, 1, 1+
  - 3: Kp = 2-, 2, 2+
  - 4: Kp = 3-, 3, 3+
  - 5: Kp = 4-, 4, 4+
  - 6: Kp = 5-, 5, 5+
  - 7: Kp >= 6-

### T96, T01, TS04
- `pdyn`: Solar wind dynamic pressure [nPa]
- `dst`: Dst index [nT]
- `byimf`: IMF By component [nT]
- `bzimf`: IMF Bz component [nT]

### Additional for T01
- `g1`, `g2`: G-indices (default: 0.0)

### Additional for TS04
- `w1` to `w6`: W-indices (default: 0.0)

## License

MIT
