"""
TsyganenkoModels Python Wrapper

A Python wrapper for the TsyganenkoModels.jl Julia package, providing
implementations of Tsyganenko empirical magnetospheric magnetic field models.

Supported models:
- T89: Tsyganenko 1989 model (Kp-dependent)
- T96: Tsyganenko 1996 model (solar wind dependent)
- T01: Tsyganenko 2001 model (with dawn-dusk asymmetry)
- TS04: Tsyganenko-Sitnov 2004 model (storm-time dynamics)

Example:
    >>> import numpy as np
    >>> from tsyganenko import t89, T89
    >>> # Single position
    >>> B = t89(-5.1, 0.3, 2.8, ps=-0.533, iopt=2)
    >>> # Array of positions
    >>> positions = np.array([[-5.1, 0.3, 2.8], [-4.0, 0.5, 3.0]])
    >>> B = t89(positions, ps=-0.533, iopt=2)
"""

import numpy as np
from datetime import datetime
from typing import Union, Tuple, List, Optional

# Lazy loading of Julia
_jl = None
_TsyganenkoModels = None


def _init_julia():
    """Initialize Julia and load TsyganenkoModels package."""
    global _jl, _TsyganenkoModels
    if _jl is None:
        from juliacall import Main as jl
        _jl = jl
        _jl.seval("using TsyganenkoModels")
        _TsyganenkoModels = _jl.TsyganenkoModels


def _ensure_julia():
    """Ensure Julia is initialized."""
    if _jl is None:
        _init_julia()


ArrayLike = Union[np.ndarray, List, Tuple]


def _to_array(positions: ArrayLike) -> np.ndarray:
    """Convert input to numpy array and validate shape.

    Args:
        positions: Position(s) in GSM coordinates. Can be:
            - 1D array/list/tuple of shape (3,) for single position
            - 2D array of shape (N, 3) for N positions

    Returns:
        2D numpy array of shape (N, 3)
    """
    arr = np.asarray(positions, dtype=np.float64)
    if arr.ndim == 1:
        if arr.shape[0] != 3:
            raise ValueError(f"1D input must have exactly 3 elements, got {arr.shape[0]}")
        arr = arr.reshape(1, 3)
    elif arr.ndim == 2:
        if arr.shape[1] != 3:
            raise ValueError(f"2D input must have shape (N, 3), got {arr.shape}")
    else:
        raise ValueError(f"Input must be 1D or 2D array, got {arr.ndim}D")
    return arr


def _call_model(model_func, positions: np.ndarray, *args, **kwargs) -> np.ndarray:
    """Call Julia model function for each position.

    Args:
        model_func: Julia model function to call
        positions: Array of shape (N, 3) with positions
        *args, **kwargs: Additional arguments for the model

    Returns:
        Array of shape (N, 3) with magnetic field components [Bx, By, Bz]
    """
    _ensure_julia()
    n_positions = positions.shape[0]
    results = np.zeros((n_positions, 3), dtype=np.float64)

    for i in range(n_positions):
        x, y, z = positions[i]
        result = model_func(float(x), float(y), float(z), *args, **kwargs)
        # Extract components from GSM result
        results[i, 0] = float(result[0])
        results[i, 1] = float(result[1])
        results[i, 2] = float(result[2])

    return results


def dipole_tilt(time: datetime) -> float:
    """Calculate the geodipole tilt angle.

    Args:
        time: Datetime for which to calculate the tilt angle

    Returns:
        Dipole tilt angle in radians
    """
    _ensure_julia()
    # Convert Python datetime to Julia DateTime
    jl_datetime = _jl.seval("using Dates; Dates.DateTime")(
        time.year, time.month, time.day,
        time.hour, time.minute, time.second
    )
    return float(_TsyganenkoModels.dipole_tilt(jl_datetime))


# =============================================================================
# Functional API
# =============================================================================

def t89(
    positions: ArrayLike,
    *,
    ps: float,
    iopt: int,
) -> np.ndarray:
    """Calculate magnetic field using Tsyganenko 1989 model.

    Args:
        positions: Position(s) in GSM coordinates [Earth radii].
            - Shape (3,) for single position [x, y, z]
            - Shape (N, 3) for N positions
        ps: Geodipole tilt angle [radians]
        iopt: Kp index level (1-7):
            - 1: Kp = 0, 0+
            - 2: Kp = 1-, 1, 1+
            - 3: Kp = 2-, 2, 2+
            - 4: Kp = 3-, 3, 3+
            - 5: Kp = 4-, 4, 4+
            - 6: Kp = 5-, 5, 5+
            - 7: Kp >= 6-

    Returns:
        Magnetic field [nT] in GSM coordinates.
            - Shape (3,) if single position input
            - Shape (N, 3) if multiple positions input

    Example:
        >>> B = t89([-5.1, 0.3, 2.8], ps=-0.533, iopt=2)
        >>> positions = np.array([[-5.1, 0.3, 2.8], [-4.0, 0.5, 3.0]])
        >>> B = t89(positions, ps=-0.533, iopt=2)
    """
    _ensure_julia()
    arr = _to_array(positions)
    single_input = np.asarray(positions).ndim == 1

    results = _call_model(_TsyganenkoModels.t89, arr, float(ps), int(iopt))

    return results[0] if single_input else results


def t96(
    positions: ArrayLike,
    *,
    ps: float,
    pdyn: float,
    dst: float,
    byimf: float,
    bzimf: float,
) -> np.ndarray:
    """Calculate magnetic field using Tsyganenko 1996 model.

    Args:
        positions: Position(s) in GSM coordinates [Earth radii].
            - Shape (3,) for single position [x, y, z]
            - Shape (N, 3) for N positions
        ps: Geodipole tilt angle [radians]
        pdyn: Solar wind dynamic pressure [nPa] (valid: 0.5-10)
        dst: Dst index [nT] (valid: -100 to +20)
        byimf: IMF By component [nT] (valid: -10 to +10)
        bzimf: IMF Bz component [nT] (valid: -10 to +10)

    Returns:
        Magnetic field [nT] in GSM coordinates.
            - Shape (3,) if single position input
            - Shape (N, 3) if multiple positions input

    Example:
        >>> B = t96([-5.1, 0.3, 2.8], ps=-0.533, pdyn=2.0, dst=-87.0, byimf=2.0, bzimf=-5.0)
    """
    _ensure_julia()
    arr = _to_array(positions)
    single_input = np.asarray(positions).ndim == 1

    results = _call_model(
        _TsyganenkoModels.t96, arr,
        float(ps), float(pdyn), float(dst), float(byimf), float(bzimf)
    )

    return results[0] if single_input else results


def t01(
    positions: ArrayLike,
    *,
    ps: float,
    pdyn: float,
    dst: float,
    byimf: float,
    bzimf: float,
    g1: float = 0.0,
    g2: float = 0.0,
) -> np.ndarray:
    """Calculate magnetic field using Tsyganenko 2001 model.

    Note: This model is based on data taken sunward from X=-15 Re,
    and becomes invalid at larger tailward distances.

    Args:
        positions: Position(s) in GSM coordinates [Earth radii].
            - Shape (3,) for single position [x, y, z]
            - Shape (N, 3) for N positions
        ps: Geodipole tilt angle [radians]
        pdyn: Solar wind dynamic pressure [nPa]
        dst: Dst index [nT]
        byimf: IMF By component [nT]
        bzimf: IMF Bz component [nT]
        g1: G1 index (default: 0.0)
        g2: G2 index (default: 0.0)

    Returns:
        Magnetic field [nT] in GSM coordinates.
            - Shape (3,) if single position input
            - Shape (N, 3) if multiple positions input

    Example:
        >>> B = t01([-5.1, 0.3, 2.8], ps=-0.533, pdyn=2.0, dst=-87.0, byimf=2.0, bzimf=-5.0)
    """
    _ensure_julia()
    arr = _to_array(positions)
    single_input = np.asarray(positions).ndim == 1

    results = _call_model(
        _TsyganenkoModels.t01, arr,
        float(ps), float(pdyn), float(dst), float(byimf), float(bzimf),
        g1=float(g1), g2=float(g2)
    )

    return results[0] if single_input else results


def ts04(
    positions: ArrayLike,
    *,
    ps: float,
    pdyn: float,
    dst: float,
    byimf: float,
    bzimf: float,
    w1: float = 0.0,
    w2: float = 0.0,
    w3: float = 0.0,
    w4: float = 0.0,
    w5: float = 0.0,
    w6: float = 0.0,
) -> np.ndarray:
    """Calculate magnetic field using Tsyganenko-Sitnov 2004 model.

    Args:
        positions: Position(s) in GSM coordinates [Earth radii].
            - Shape (3,) for single position [x, y, z]
            - Shape (N, 3) for N positions
        ps: Geodipole tilt angle [radians]
        pdyn: Solar wind dynamic pressure [nPa]
        dst: Dst index [nT]
        byimf: IMF By component [nT]
        bzimf: IMF Bz component [nT]
        w1-w6: W indices (default: 0.0)

    Returns:
        Magnetic field [nT] in GSM coordinates.
            - Shape (3,) if single position input
            - Shape (N, 3) if multiple positions input

    Example:
        >>> B = ts04([-5.1, 0.3, 2.8], ps=-0.533, pdyn=2.0, dst=-87.0, byimf=2.0, bzimf=-5.0)
    """
    _ensure_julia()
    arr = _to_array(positions)
    single_input = np.asarray(positions).ndim == 1

    results = _call_model(
        _TsyganenkoModels.ts04, arr,
        float(ps), float(pdyn), float(dst), float(byimf), float(bzimf),
        w1=float(w1), w2=float(w2), w3=float(w3),
        w4=float(w4), w5=float(w5), w6=float(w6)
    )

    return results[0] if single_input else results


# =============================================================================
# Object-Oriented API
# =============================================================================

class T89:
    """Tsyganenko 1989 model.

    Args:
        iopt: Kp index level (1-7)

    Example:
        >>> model = T89(iopt=2)
        >>> B = model([-5.1, 0.3, 2.8], ps=-0.533)
    """

    def __init__(self, iopt: int):
        if not 1 <= iopt <= 7:
            raise ValueError("iopt must be between 1 and 7")
        self.iopt = iopt
        self._jl_model = None

    def _get_model(self):
        """Get or create Julia model instance."""
        _ensure_julia()
        if self._jl_model is None:
            self._jl_model = _TsyganenkoModels.T89(self.iopt)
        return self._jl_model

    def __call__(
        self,
        positions: ArrayLike,
        ps: float,
    ) -> np.ndarray:
        """Compute magnetic field at given position(s).

        Args:
            positions: Position(s) in GSM coordinates [Earth radii]
            ps: Geodipole tilt angle [radians]

        Returns:
            Magnetic field [nT] in GSM coordinates
        """
        arr = _to_array(positions)
        single_input = np.asarray(positions).ndim == 1
        model = self._get_model()

        n_positions = arr.shape[0]
        results = np.zeros((n_positions, 3), dtype=np.float64)

        for i in range(n_positions):
            pos = tuple(arr[i])
            result = model(pos, float(ps))
            results[i, 0] = float(result[0])
            results[i, 1] = float(result[1])
            results[i, 2] = float(result[2])

        return results[0] if single_input else results

    def __repr__(self):
        return f"T89(iopt={self.iopt})"


class T96:
    """Tsyganenko 1996 model.

    Args:
        pdyn: Solar wind dynamic pressure [nPa]
        dst: Dst index [nT]
        byimf: IMF By component [nT]
        bzimf: IMF Bz component [nT]

    Example:
        >>> model = T96(pdyn=2.0, dst=-87.0, byimf=2.0, bzimf=-5.0)
        >>> B = model([-5.1, 0.3, 2.8], ps=-0.533)
    """

    def __init__(self, pdyn: float, dst: float, byimf: float, bzimf: float):
        self.pdyn = pdyn
        self.dst = dst
        self.byimf = byimf
        self.bzimf = bzimf
        self._jl_model = None

    def _get_model(self):
        """Get or create Julia model instance."""
        _ensure_julia()
        if self._jl_model is None:
            self._jl_model = _TsyganenkoModels.T96(
                pdyn=self.pdyn, dst=self.dst,
                byimf=self.byimf, bzimf=self.bzimf
            )
        return self._jl_model

    def __call__(
        self,
        positions: ArrayLike,
        ps: float,
    ) -> np.ndarray:
        """Compute magnetic field at given position(s).

        Args:
            positions: Position(s) in GSM coordinates [Earth radii]
            ps: Geodipole tilt angle [radians]

        Returns:
            Magnetic field [nT] in GSM coordinates
        """
        arr = _to_array(positions)
        single_input = np.asarray(positions).ndim == 1
        model = self._get_model()

        n_positions = arr.shape[0]
        results = np.zeros((n_positions, 3), dtype=np.float64)

        for i in range(n_positions):
            pos = tuple(arr[i])
            result = model(pos, float(ps))
            results[i, 0] = float(result[0])
            results[i, 1] = float(result[1])
            results[i, 2] = float(result[2])

        return results[0] if single_input else results

    def __repr__(self):
        return f"T96(pdyn={self.pdyn}, dst={self.dst}, byimf={self.byimf}, bzimf={self.bzimf})"


class T01:
    """Tsyganenko 2001 model.

    Note: This model is based on data taken sunward from X=-15 Re,
    and becomes invalid at larger tailward distances.

    Args:
        pdyn: Solar wind dynamic pressure [nPa]
        dst: Dst index [nT]
        byimf: IMF By component [nT]
        bzimf: IMF Bz component [nT]
        g1: G1 index (default: 0.0)
        g2: G2 index (default: 0.0)

    Example:
        >>> model = T01(pdyn=2.0, dst=-87.0, byimf=2.0, bzimf=-5.0)
        >>> B = model([-5.1, 0.3, 2.8], ps=-0.533)
    """

    def __init__(
        self,
        pdyn: float,
        dst: float,
        byimf: float,
        bzimf: float,
        g1: float = 0.0,
        g2: float = 0.0,
    ):
        self.pdyn = pdyn
        self.dst = dst
        self.byimf = byimf
        self.bzimf = bzimf
        self.g1 = g1
        self.g2 = g2
        self._jl_model = None

    def _get_model(self):
        """Get or create Julia model instance."""
        _ensure_julia()
        if self._jl_model is None:
            self._jl_model = _TsyganenkoModels.T01(
                pdyn=self.pdyn, dst=self.dst,
                byimf=self.byimf, bzimf=self.bzimf,
                g1=self.g1, g2=self.g2
            )
        return self._jl_model

    def __call__(
        self,
        positions: ArrayLike,
        ps: float,
    ) -> np.ndarray:
        """Compute magnetic field at given position(s).

        Args:
            positions: Position(s) in GSM coordinates [Earth radii]
            ps: Geodipole tilt angle [radians]

        Returns:
            Magnetic field [nT] in GSM coordinates
        """
        arr = _to_array(positions)
        single_input = np.asarray(positions).ndim == 1
        model = self._get_model()

        n_positions = arr.shape[0]
        results = np.zeros((n_positions, 3), dtype=np.float64)

        for i in range(n_positions):
            pos = tuple(arr[i])
            result = model(pos, float(ps))
            results[i, 0] = float(result[0])
            results[i, 1] = float(result[1])
            results[i, 2] = float(result[2])

        return results[0] if single_input else results

    def __repr__(self):
        return f"T01(pdyn={self.pdyn}, dst={self.dst}, byimf={self.byimf}, bzimf={self.bzimf}, g1={self.g1}, g2={self.g2})"


class TS04:
    """Tsyganenko-Sitnov 2004 model.

    Args:
        pdyn: Solar wind dynamic pressure [nPa]
        dst: Dst index [nT]
        byimf: IMF By component [nT]
        bzimf: IMF Bz component [nT]
        w1-w6: W indices (default: 0.0)

    Example:
        >>> model = TS04(pdyn=2.0, dst=-87.0, byimf=2.0, bzimf=-5.0)
        >>> B = model([-5.1, 0.3, 2.8], ps=-0.533)
    """

    def __init__(
        self,
        pdyn: float,
        dst: float,
        byimf: float,
        bzimf: float,
        w1: float = 0.0,
        w2: float = 0.0,
        w3: float = 0.0,
        w4: float = 0.0,
        w5: float = 0.0,
        w6: float = 0.0,
    ):
        self.pdyn = pdyn
        self.dst = dst
        self.byimf = byimf
        self.bzimf = bzimf
        self.w1 = w1
        self.w2 = w2
        self.w3 = w3
        self.w4 = w4
        self.w5 = w5
        self.w6 = w6
        self._jl_model = None

    def _get_model(self):
        """Get or create Julia model instance."""
        _ensure_julia()
        if self._jl_model is None:
            self._jl_model = _TsyganenkoModels.TS04(
                pdyn=self.pdyn, dst=self.dst,
                byimf=self.byimf, bzimf=self.bzimf,
                w1=self.w1, w2=self.w2, w3=self.w3,
                w4=self.w4, w5=self.w5, w6=self.w6
            )
        return self._jl_model

    def __call__(
        self,
        positions: ArrayLike,
        ps: float,
    ) -> np.ndarray:
        """Compute magnetic field at given position(s).

        Args:
            positions: Position(s) in GSM coordinates [Earth radii]
            ps: Geodipole tilt angle [radians]

        Returns:
            Magnetic field [nT] in GSM coordinates
        """
        arr = _to_array(positions)
        single_input = np.asarray(positions).ndim == 1
        model = self._get_model()

        n_positions = arr.shape[0]
        results = np.zeros((n_positions, 3), dtype=np.float64)

        for i in range(n_positions):
            pos = tuple(arr[i])
            result = model(pos, float(ps))
            results[i, 0] = float(result[0])
            results[i, 1] = float(result[1])
            results[i, 2] = float(result[2])

        return results[0] if single_input else results

    def __repr__(self):
        return f"TS04(pdyn={self.pdyn}, dst={self.dst}, byimf={self.byimf}, bzimf={self.bzimf}, w1={self.w1}, w2={self.w2}, w3={self.w3}, w4={self.w4}, w5={self.w5}, w6={self.w6})"


__all__ = [
    # Functions
    "t89", "t96", "t01", "ts04",
    "dipole_tilt",
    # Classes
    "T89", "T96", "T01", "TS04",
]

__version__ = "0.1.0"
