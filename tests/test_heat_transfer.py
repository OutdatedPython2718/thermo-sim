"""Tests for heat transfer simulation modules."""

import numpy as np
import pytest
from scipy.special import erfc

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


class TestTransientSemiInfiniteSolid:
    """Benchmark A: Semi-infinite solid with sudden surface temperature.

    Analytical: T(x,t) = T_surface * erfc(x / (2 * sqrt(alpha * t)))
    Ref: Incropera & DeWitt, Ch. 5, Eq. 5.60
    """

    def test_erfc_benchmark(self):
        T_surface = 100.0
        T_initial_val = 0.0
        alpha = 1e-5
        L = 0.5  # long enough that far boundary doesn't matter

        nx = 201
        dx = L / (nx - 1)
        dt = 0.4 * dx**2 / alpha  # Fo = 0.4
        t_final = 50.0
        n_steps = int(t_final / dt)

        T_init = np.full(nx, T_initial_val)
        history = solve_1d_transient_explicit(
            T_initial=T_init, L=L, alpha=alpha,
            dt=dt, n_steps=n_steps,
            T_left=T_surface, T_right=T_initial_val,
            save_every=n_steps,
        )
        T_numerical = history[-1]

        x = np.linspace(0, L, nx)
        T_analytical = T_surface * erfc(x / (2 * np.sqrt(alpha * t_final)))

        check_indices = np.linspace(1, nx - 2, 5, dtype=int)
        for idx in check_indices:
            if T_analytical[idx] < 1.0:
                continue
            rel_err = abs(T_numerical[idx] - T_analytical[idx]) / T_analytical[idx]
            assert rel_err < 0.02, (
                f"Semi-infinite solid at x={x[idx]:.3f}m: "
                f"numerical={T_numerical[idx]:.4f}, "
                f"analytical={T_analytical[idx]:.4f}, "
                f"rel_err={rel_err:.4f}"
            )


class TestTransientFourierSeries:
    """Benchmark B: Finite bar with sudden boundary temperature change.

    Analytical: T(x,t) = T_hot*(1 - x/L) - sum (2*T_hot/(n*pi))
                * sin(n*pi*x/L) * exp(-n^2*pi^2*alpha*t/L^2)
    Ref: Carslaw & Jaeger, Section 3.3
    """

    @staticmethod
    def analytical_fourier(x, t, L, alpha, T_hot, n_terms=50):
        T = T_hot * (1.0 - x / L)
        for n in range(1, n_terms + 1):
            T -= (
                (2 * T_hot / (n * np.pi))
                * np.sin(n * np.pi * x / L)
                * np.exp(-n**2 * np.pi**2 * alpha * t / L**2)
            )
        return T

    def test_moderate_time(self):
        T_hot = 100.0
        L = 1.0
        alpha = 1e-4
        nx = 101
        dx = L / (nx - 1)
        dt = 0.4 * dx**2 / alpha
        t_final = 200.0
        n_steps = int(t_final / dt)

        T_init = np.zeros(nx)
        history = solve_1d_transient_explicit(
            T_initial=T_init, L=L, alpha=alpha,
            dt=dt, n_steps=n_steps,
            T_left=T_hot, T_right=0.0,
            save_every=n_steps,
        )
        T_numerical = history[-1]

        x = np.linspace(0, L, nx)
        T_analytical = self.analytical_fourier(x, t_final, L, alpha, T_hot)

        check_indices = np.linspace(10, nx - 10, 5, dtype=int)
        for idx in check_indices:
            if T_analytical[idx] < 1.0:
                continue
            rel_err = abs(T_numerical[idx] - T_analytical[idx]) / T_analytical[idx]
            assert rel_err < 0.02, (
                f"Fourier bar at x={x[idx]:.3f}m, t={t_final}s: "
                f"numerical={T_numerical[idx]:.4f}, "
                f"analytical={T_analytical[idx]:.4f}, "
                f"rel_err={rel_err:.4f}"
            )

    def test_long_time_matches_steady_state(self):
        T_hot = 100.0
        L = 1.0
        alpha = 1e-4
        nx = 51
        dx = L / (nx - 1)
        dt = 0.4 * dx**2 / alpha
        t_final = 5000.0
        n_steps = int(t_final / dt)

        T_init = np.zeros(nx)
        history = solve_1d_transient_explicit(
            T_initial=T_init, L=L, alpha=alpha,
            dt=dt, n_steps=n_steps,
            T_left=T_hot, T_right=0.0,
            save_every=n_steps,
        )
        T_numerical = history[-1]

        x = np.linspace(0, L, nx)
        T_steady = T_hot * (1.0 - x / L)

        for idx in range(1, nx - 1):
            if T_steady[idx] < 1.0:
                continue
            rel_err = abs(T_numerical[idx] - T_steady[idx]) / T_steady[idx]
            assert rel_err < 0.02, (
                f"Steady-state at x={x[idx]:.3f}m: "
                f"numerical={T_numerical[idx]:.4f}, "
                f"steady={T_steady[idx]:.4f}, "
                f"rel_err={rel_err:.4f}"
            )


class TestGridConvergenceOrder:
    """Verify 2D Gauss-Seidel solver converges at 2nd order.

    Probes at the quarter-point (n//4, n//4) to avoid the
    domain center, where symmetric BCs produce exact values
    regardless of grid resolution.
    """

    def test_observed_order_is_second(self):
        bc = BoundaryCondition(
            top=400.0, bottom=200.0, left=300.0, right=100.0,
        )
        resolutions = [21, 41, 81]
        probe_temps = []
        dx_values = []

        for n in resolutions:
            T = solve_2d_steady(
                nx=n, ny=n, Lx=1.0, Ly=1.0,
                bc=bc, k=1.0,
                tol=1e-10, max_iter=200000,
            )
            # Probe at quarter-point to avoid center symmetry
            j = n // 4
            i = n // 4
            probe_temps.append(T[j, i])
            dx_values.append(1.0 / (n - 1))

        # Richardson extrapolation for reference
        T_exact = (
            probe_temps[-1]
            + (probe_temps[-1] - probe_temps[-2]) / 3
        )
        errors = [abs(tc - T_exact) for tc in probe_temps]

        for i in range(len(errors) - 1):
            if errors[i + 1] < 1e-12:
                continue
            p = np.log(errors[i] / errors[i + 1]) / np.log(
                dx_values[i] / dx_values[i + 1]
            )
            assert 1.8 < p < 2.2, (
                f"Convergence order between "
                f"{resolutions[i]}x{resolutions[i]} and "
                f"{resolutions[i+1]}x{resolutions[i+1]}: "
                f"p={p:.2f}, expected ~2.0"
            )
