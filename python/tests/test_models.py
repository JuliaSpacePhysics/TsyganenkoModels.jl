"""Tests for TsyganenkoModels Python wrapper."""

import numpy as np
import pytest
from datetime import datetime

from tsyganenko import t89, t96, t01, ts04, dipole_tilt
from tsyganenko import T89, T96, T01, TS04


# Reference values from Julia tests
PS = -0.533585131
POSITION = (-5.1, 0.3, 2.8)

# Expected values from TsyganenkoModels.jl test suite
T89_EXPECTED = np.array([20.774513266271535, -0.6468620203393424, -15.076847577020597])
T96_EXPECTED = np.array([61.17527804691476, -1.4613918028982884, -40.44600024617623])


class TestArrayInput:
    """Test array input handling."""

    def test_list_input(self):
        """Test list input."""
        B = t89([-5.1, 0.3, 2.8], ps=PS, iopt=2)
        assert B.shape == (3,)

    def test_tuple_input(self):
        """Test tuple input."""
        B = t89((-5.1, 0.3, 2.8), ps=PS, iopt=2)
        assert B.shape == (3,)

    def test_numpy_1d_input(self):
        """Test 1D numpy array input."""
        B = t89(np.array([-5.1, 0.3, 2.8]), ps=PS, iopt=2)
        assert B.shape == (3,)

    def test_numpy_2d_input(self):
        """Test 2D numpy array input."""
        positions = np.array([
            [-5.1, 0.3, 2.8],
            [-4.0, 0.5, 3.0],
        ])
        B = t89(positions, ps=PS, iopt=2)
        assert B.shape == (2, 3)

    def test_list_of_lists_input(self):
        """Test list of lists input."""
        positions = [[-5.1, 0.3, 2.8], [-4.0, 0.5, 3.0]]
        B = t89(positions, ps=PS, iopt=2)
        assert B.shape == (2, 3)

    def test_invalid_shape_raises(self):
        """Test that invalid shapes raise ValueError."""
        with pytest.raises(ValueError):
            t89([1, 2], ps=PS, iopt=2)  # Wrong length

        with pytest.raises(ValueError):
            t89(np.array([[1, 2], [3, 4]]), ps=PS, iopt=2)  # Wrong shape


class TestT89:
    """Test T89 model."""

    def test_functional_api(self):
        """Test functional API returns correct values."""
        B = t89(POSITION, ps=PS, iopt=2)
        np.testing.assert_allclose(B, T89_EXPECTED, rtol=1e-10)

    def test_class_api(self):
        """Test class API returns correct values."""
        model = T89(iopt=2)
        B = model(POSITION, ps=PS)
        np.testing.assert_allclose(B, T89_EXPECTED, rtol=1e-10)

    def test_array_input(self):
        """Test array input returns correct shape."""
        positions = np.array([POSITION, POSITION])
        B = t89(positions, ps=PS, iopt=2)
        assert B.shape == (2, 3)
        np.testing.assert_allclose(B[0], T89_EXPECTED, rtol=1e-10)
        np.testing.assert_allclose(B[1], T89_EXPECTED, rtol=1e-10)

    def test_iopt_range(self):
        """Test iopt validation."""
        with pytest.raises(ValueError):
            T89(iopt=0)
        with pytest.raises(ValueError):
            T89(iopt=8)

    def test_repr(self):
        """Test model representation."""
        model = T89(iopt=2)
        assert "T89(iopt=2)" == repr(model)


class TestT96:
    """Test T96 model."""

    def test_functional_api(self):
        """Test functional API returns correct values."""
        B = t96(POSITION, ps=PS, pdyn=2.0, dst=-87.0, byimf=2.0, bzimf=-5.0)
        np.testing.assert_allclose(B, T96_EXPECTED, rtol=1e-10)

    def test_class_api(self):
        """Test class API returns correct values."""
        model = T96(pdyn=2.0, dst=-87.0, byimf=2.0, bzimf=-5.0)
        B = model(POSITION, ps=PS)
        np.testing.assert_allclose(B, T96_EXPECTED, rtol=1e-10)

    def test_array_input(self):
        """Test array input returns correct shape."""
        positions = np.array([POSITION, POSITION])
        model = T96(pdyn=2.0, dst=-87.0, byimf=2.0, bzimf=-5.0)
        B = model(positions, ps=PS)
        assert B.shape == (2, 3)

    def test_repr(self):
        """Test model representation."""
        model = T96(pdyn=2.0, dst=-87.0, byimf=2.0, bzimf=-5.0)
        assert "pdyn=2.0" in repr(model)


class TestT01:
    """Test T01 model."""

    def test_functional_api(self):
        """Test functional API."""
        B = t01(POSITION, ps=PS, pdyn=2.0, dst=-87.0, byimf=2.0, bzimf=-5.0)
        assert B.shape == (3,)

    def test_class_api(self):
        """Test class API."""
        model = T01(pdyn=2.0, dst=-87.0, byimf=2.0, bzimf=-5.0)
        B = model(POSITION, ps=PS)
        assert B.shape == (3,)

    def test_g_indices(self):
        """Test G-indices parameters."""
        model = T01(pdyn=2.0, dst=-87.0, byimf=2.0, bzimf=-5.0, g1=1.0, g2=0.5)
        B = model(POSITION, ps=PS)
        assert B.shape == (3,)

    def test_repr(self):
        """Test model representation."""
        model = T01(pdyn=2.0, dst=-87.0, byimf=2.0, bzimf=-5.0, g1=1.0, g2=0.5)
        assert "g1=1.0" in repr(model)
        assert "g2=0.5" in repr(model)


class TestTS04:
    """Test TS04 model."""

    def test_functional_api(self):
        """Test functional API."""
        B = ts04(POSITION, ps=PS, pdyn=2.0, dst=-87.0, byimf=2.0, bzimf=-5.0)
        assert B.shape == (3,)

    def test_class_api(self):
        """Test class API."""
        model = TS04(pdyn=2.0, dst=-87.0, byimf=2.0, bzimf=-5.0)
        B = model(POSITION, ps=PS)
        assert B.shape == (3,)

    def test_w_indices(self):
        """Test W-indices parameters."""
        model = TS04(
            pdyn=2.0, dst=-87.0, byimf=2.0, bzimf=-5.0,
            w1=0.5, w2=0.3, w3=0.2, w4=0.1, w5=0.05, w6=0.01
        )
        B = model(POSITION, ps=PS)
        assert B.shape == (3,)

    def test_repr(self):
        """Test model representation."""
        model = TS04(pdyn=2.0, dst=-87.0, byimf=2.0, bzimf=-5.0, w1=0.5)
        assert "w1=0.5" in repr(model)


class TestDipoleTilt:
    """Test dipole tilt calculation."""

    def test_dipole_tilt(self):
        """Test dipole tilt calculation returns reasonable value."""
        time = datetime(2001, 1, 1, 2, 3, 4)
        ps = dipole_tilt(time)
        # Dipole tilt should be in reasonable range (-0.6 to 0.6 radians)
        assert -0.6 < ps < 0.6


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
