"""Tests for heat transfer simulation modules."""

import numpy as np
import pytest

from simulations.heat_transfer.conduction import (
    BoundaryCondition,
    create_grid,
    solve_1d_steady,
    solve_2d_steady,
    solve_2d_transient_explicit,
)
from simulations.heat_transfer.conduction import solve_1d_transient_explicit


class TestGrid:
    def test_create_grid(self):
        grid = create_grid(Lx=1.0, Ly=0.5, nx=11, ny=6)
        assert grid["x"].shape == (11,)
        assert grid["y"].shape == (6,)
        assert grid["dx"] == pytest.approx(0.1)
        assert grid["dy"] == pytest.approx(0.1)


class TestSteady1D:
    def test_linear_profile(self):
        T = solve_1d_steady(nx=51, L=1.0, T_left=100.0, T_right=200.0, k=1.0)
        assert len(T) == 51
        assert T[0] == pytest.approx(100.0)
        assert T[-1] == pytest.approx(200.0)
        assert T[25] == pytest.approx(150.0, abs=1.0)

    def test_with_source(self):
        T_no_source = solve_1d_steady(nx=51, L=1.0, T_left=100.0, T_right=100.0, k=1.0)
        source = np.ones(51) * 1000.0
        T_source = solve_1d_steady(nx=51, L=1.0, T_left=100.0, T_right=100.0, k=1.0, source=source)
        assert T_source[25] > T_no_source[25]


class TestSteady2D:
    def test_uniform_boundary(self):
        bc = BoundaryCondition(top=100.0, bottom=100.0, left=100.0, right=100.0)
        T = solve_2d_steady(nx=11, ny=11, Lx=1.0, Ly=1.0, bc=bc, k=1.0)
        assert T.shape == (11, 11)
        assert np.allclose(T, 100.0, atol=0.1)

    def test_gradient(self):
        bc = BoundaryCondition(top=200.0, bottom=0.0, left=None, right=None)
        T = solve_2d_steady(nx=11, ny=21, Lx=1.0, Ly=1.0, bc=bc, k=1.0)
        assert T[10, 5] == pytest.approx(100.0, abs=10.0)
        assert T[0, 5] > T[10, 5] > T[-1, 5]


class TestTransient:
    def test_approaches_steady_state(self):
        bc = BoundaryCondition(top=100.0, bottom=0.0, left=0.0, right=0.0)
        T_initial = np.zeros((21, 21))
        T_initial[0, :] = 100.0
        T_history = solve_2d_transient_explicit(
            T_initial=T_initial, Lx=1.0, Ly=1.0, alpha=1e-4,
            dt=0.1, n_steps=5000, bc=bc,
        )
        T_final = T_history[-1]
        assert 0 < T_final[10, 10] < 100

    def test_fourier_number_warning(self):
        bc = BoundaryCondition(top=100.0, bottom=0.0, left=0.0, right=0.0)
        T_initial = np.zeros((11, 11))
        T_initial[0, :] = 100.0
        with pytest.warns(UserWarning, match="Fourier number"):
            solve_2d_transient_explicit(
                T_initial=T_initial, Lx=1.0, Ly=1.0, alpha=1.0,
                dt=1.0, n_steps=1, bc=bc,
            )


class TestTransient1D:
    def test_returns_history(self):
        nx = 21
        T_initial = np.zeros(nx)
        T_initial[0] = 100.0
        history = solve_1d_transient_explicit(
            T_initial=T_initial, L=1.0, alpha=1e-4,
            dt=0.1, n_steps=100, T_left=100.0, T_right=0.0,
            save_every=50,
        )
        assert len(history) >= 2
        assert len(history[0]) == nx
        assert history[0][0] == pytest.approx(100.0)
