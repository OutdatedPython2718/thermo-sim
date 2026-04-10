"""Numerical solvers for thermodynamic simulations.

Provides RK4 integration, fixed-point iterative solving,
and Gauss-Seidel 2D steady-state heat conduction.
"""

from dataclasses import dataclass

import numpy as np


def rk4_integrate(f, y0, t_span, dt):
    """4th-order Runge-Kutta integration.

    Parameters
    ----------
    f : callable — f(t, y) -> dy/dt
    y0 : float or ndarray — initial condition
    t_span : tuple — (t_start, t_end)
    dt : float — time step

    Returns
    -------
    t_arr, y_arr : ndarrays
    """
    t_start, t_end = t_span
    n_steps = int(round((t_end - t_start) / dt)) + 1
    t_arr = np.linspace(t_start, t_end, n_steps)

    scalar = np.isscalar(y0)
    y = np.float64(y0) if scalar else np.array(y0, dtype=np.float64)

    if scalar:
        y_arr = np.zeros(n_steps)
    else:
        y_arr = np.zeros((n_steps, len(y)))

    y_arr[0] = y

    for i in range(n_steps - 1):
        t = t_arr[i]
        h = t_arr[i + 1] - t_arr[i]
        k1 = np.asarray(f(t, y), dtype=np.float64)
        k2 = np.asarray(f(t + h / 2, y + h / 2 * k1), dtype=np.float64)
        k3 = np.asarray(f(t + h / 2, y + h / 2 * k2), dtype=np.float64)
        k4 = np.asarray(f(t + h, y + h * k3), dtype=np.float64)
        y = y + (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
        y_arr[i + 1] = y

    return t_arr, y_arr


@dataclass
class IterativeResult:
    """Result from iterative_solve."""
    x: float | np.ndarray
    converged: bool
    iterations: int
    residual: float


def iterative_solve(g, x0, tol=1e-8, max_iter=200):
    """Fixed-point iteration: find x such that x = g(x)."""
    x = np.float64(x0) if np.isscalar(x0) else np.array(x0, dtype=np.float64)

    for i in range(max_iter):
        x_new = g(x)
        residual = float(np.max(np.abs(x_new - x)))
        x = x_new if np.isscalar(x_new) else np.array(x_new, dtype=np.float64)
        if residual < tol:
            return IterativeResult(x=x, converged=True, iterations=i + 1, residual=residual)
        if not np.all(np.isfinite(x)):
            return IterativeResult(x=x, converged=False, iterations=i + 1, residual=float("inf"))

    return IterativeResult(x=x, converged=False, iterations=max_iter, residual=residual)


def gauss_seidel_2d(
    nx, ny, bc_top=0.0, bc_bottom=0.0, bc_left=0.0, bc_right=0.0,
    tol=1e-6, max_iter=10000, source=None
):
    """Gauss-Seidel solver for 2D steady-state heat conduction (Laplace/Poisson)."""
    T = np.zeros((ny, nx))
    T[0, :] = bc_top
    T[-1, :] = bc_bottom
    T[:, 0] = bc_left
    T[:, -1] = bc_right

    if source is None:
        source = np.zeros((ny, nx))

    dx = 1.0 / (nx - 1)

    for iteration in range(max_iter):
        max_change = 0.0
        for j in range(1, ny - 1):
            for i in range(1, nx - 1):
                T_old = T[j, i]
                T[j, i] = 0.25 * (
                    T[j + 1, i] + T[j - 1, i] + T[j, i + 1] + T[j, i - 1] + dx**2 * source[j, i]
                )
                change = abs(T[j, i] - T_old)
                if change > max_change:
                    max_change = change
        if max_change < tol:
            break

    return T
