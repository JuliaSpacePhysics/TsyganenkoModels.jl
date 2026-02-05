"""Python wrapper for TsyganenkoModels.jl."""

from juliacall import Main as jl

# Load the Julia package
jl.seval("using TsyganenkoModels")

# Re-export Julia functions
t89 = jl.TsyganenkoModels.t89
t96 = jl.TsyganenkoModels.t96
t01 = jl.TsyganenkoModels.t01
ts04 = jl.TsyganenkoModels.ts04
dipole_tilt = jl.TsyganenkoModels.dipole_tilt

# Model constructors/types
T89 = jl.TsyganenkoModels.T89
T96 = jl.TsyganenkoModels.T96
T01 = jl.TsyganenkoModels.T01
TS04 = jl.TsyganenkoModels.TS04
TsyIGRF = jl.TsyganenkoModels.TsyIGRF
IGRF = jl.TsyganenkoModels.IGRF

__version__ = "0.1.0"

__all__ = [
    "t89",
    "t96",
    "t01",
    "ts04",
    "dipole_tilt",
    "T89",
    "T96",
    "T01",
    "TS04",
    "TsyIGRF",
    "IGRF",
    "__version__",
    "jl",
]
