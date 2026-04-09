"""Tests for thermosim.solvers — numerical methods."""

import numpy as np
import pytest
from thermosim.solvers import rk4_integrate, iterative_solve, gauss_seidel_2d


class TestRK4:
    def test_exponential_decay(self):
        """dy/dt = -y, y(0) = 1 => y(t) = exp(-t)."""
        def f(t, y):
            return -y
        t, y = rk4_integrate(f, y0=1.0, t_span=(0.0, 2.0), dt=0.01)
        assert y[-1] == pytest.approx(np.exp(-2.0), rel=1e-6)
        assert t[-1] == pytest.approx(2.0, abs=0.01)

    def test_simple_harmonic(self):
        """dy1/dt = y2, dy2/dt = -y1 => sin/cos oscillation."""
        def f(t, y):
            return np.array([y[1], -y[0]])
        t, y = rk4_integrate(f, y0=np.array([0.0, 1.0]), t_span=(0.0, 2 * np.pi), dt=0.01)
        assert y[-1][0] == pytest.approx(0.0, abs=1e-4)
        assert y[-1][1] == pytest.approx(1.0, abs=1e-4)

    def test_returns_all_steps(self):
        def f(t, y):
            return -y
        t, y = rk4_integrate(f, y0=1.0, t_span=(0.0, 1.0), dt=0.1)
        assert len(t) == len(y)
        assert len(t) == pytest.approx(11, abs=1)


class TestIterativeSolve:
    def test_linear_system_converges(self):
        """Solve x = cos(x) iteratively. Solution ~ 0.7391."""
        def g(x):
            return np.cos(x)
        result = iterative_solve(g, x0=0.5, tol=1e-10, max_iter=100)
        assert result.converged is True
        assert result.x == pytest.approx(0.7390851332, rel=1e-8)
        assert result.iterations < 100

    def test_divergent_raises(self):
        def g(x):
            return 2 * x + 1
        result = iterative_solve(g, x0=1.0, tol=1e-10, max_iter=50)
        assert result.converged is False


class TestGaussSeidel2D:
    def test_uniform_boundary(self):
        T = gauss_seidel_2d(nx=10, ny=10, bc_top=100.0, bc_bottom=100.0, bc_left=100.0, bc_right=100.0, tol=1e-6, max_iter=1000)
        assert T.shape == (10, 10)
        assert np.allclose(T, 100.0, atol=0.1)

    def test_one_hot_boundary(self):
        T = gauss_seidel_2d(nx=20, ny=20, bc_top=100.0, bc_bottom=0.0, bc_left=0.0, bc_right=0.0, tol=1e-6, max_iter=5000)
        center = T[10, 10]
        assert 0 < center < 100
        assert T[0, 10] == pytest.approx(100.0, abs=5.0)
        assert T[-1, 10] == pytest.approx(0.0, abs=5.0)
