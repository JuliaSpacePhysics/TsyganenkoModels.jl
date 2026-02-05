import pytest


def test_import():
    import tsyganenkomodels as tm

    assert hasattr(tm, "T96")
    assert hasattr(tm, "t96")
    assert hasattr(tm, "dipole_tilt")


def test_basic_eval():
    import tsyganenkomodels as tm
    import datetime

    ps = tm.dipole_tilt(datetime.datetime(2001, 1, 1, 2, 3, 4))
    model = tm.T96(pdyn=2.0, dst=-87.0, byimf=2.0, bzimf=-5.0)
    b = model([-5.1,0.3,2.8], ps)
    assert len(b) == 3
