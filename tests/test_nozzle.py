"""Tests for nozzle flow solver."""

import numpy as np
import pytest

from simulations.nozzle_flow.nozzle import (
    NozzleGeometry,
    isentropic_mach_from_area_ratio,
    solve_nozzle_flow,
)


class TestIsentropicRelations:
    def test_sonic_at_throat(self):
        M = isentropic_mach_from_area_ratio(1.0, gamma=1.4, supersonic=False)
        assert M == pytest.approx(1.0, abs=0.01)

    def test_subsonic_for_small_area_ratio(self):
        M = isentropic_mach_from_area_ratio(2.0, gamma=1.4, supersonic=False)
        assert 0 < M < 1

    def test_supersonic_for_large_area_ratio(self):
        M = isentropic_mach_from_area_ratio(2.0, gamma=1.4, supersonic=True)
        assert M > 1


class TestNozzleFlow:
    def test_solve_produces_arrays(self):
        geom = NozzleGeometry(L=1.0, A_inlet=0.2, A_throat=0.1, A_exit=0.3)
        result = solve_nozzle_flow(geom=geom, P0=500e3, T0=300.0, gamma=1.4, R=287.0, n_points=100)
        assert len(result["x"]) == 100
        assert len(result["Mach"]) == 100

    def test_mach_one_near_throat(self):
        geom = NozzleGeometry(L=1.0, A_inlet=0.2, A_throat=0.1, A_exit=0.3)
        result = solve_nozzle_flow(geom=geom, P0=500e3, T0=300.0, gamma=1.4, R=287.0, n_points=200)
        throat_idx = np.argmin(result["A"])
        assert result["Mach"][throat_idx] == pytest.approx(1.0, abs=0.05)

    def test_pressure_decreases_in_diverging(self):
        geom = NozzleGeometry(L=1.0, A_inlet=0.2, A_throat=0.1, A_exit=0.3)
        result = solve_nozzle_flow(geom=geom, P0=500e3, T0=300.0, gamma=1.4, R=287.0, n_points=200)
        throat_idx = np.argmin(result["A"])
        assert result["P"][-1] < result["P"][throat_idx]
