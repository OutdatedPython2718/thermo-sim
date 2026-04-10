"""Finite-difference conduction solvers for 1D and 2D heat transfer."""

from __future__ import annotations

import warnings
from dataclasses import dataclass

import numpy as np


@dataclass
class BoundaryCondition:
    """Dirichlet boundary conditions. None means insulated (zero-flux)."""
    top: float | None = 0.0
    bottom: float | None = 0.0
    left: float | None = 0.0
    right: float | None = 0.0


def create_grid(Lx, Ly, nx, ny):
    dx = Lx / (nx - 1)
    dy = Ly / (ny - 1)
    x = np.linspace(0, Lx, nx)
    y = np.linspace(0, Ly, ny)
    return {"x": x, "y": y, "dx": dx, "dy": dy}


def solve_1d_steady(nx, L, T_left, T_right, k=1.0, source=None):
    """Solve 1D steady-state conduction: d/dx(k dT/dx) + source = 0."""
    dx = L / (nx - 1)
    if source is None:
        source = np.zeros(nx)
    n = nx - 2
    diag = np.full(n, -2.0)
    off_diag = np.ones(n - 1)
    rhs = -source[1:-1] * dx**2 / k
    rhs[0] -= T_left
    rhs[-1] -= T_right
    A = np.diag(diag) + np.diag(off_diag, 1) + np.diag(off_diag, -1)
    T_interior = np.linalg.solve(A, rhs)
    T = np.zeros(nx)
    T[0] = T_left
    T[-1] = T_right
    T[1:-1] = T_interior
    return T


def solve_2d_steady(nx, ny, Lx, Ly, bc, k=1.0, source=None, tol=1e-6, max_iter=10000):
    """Solve 2D steady-state conduction via Gauss-Seidel iteration."""
    dx = Lx / (nx - 1)
    dy = Ly / (ny - 1)
    T = np.zeros((ny, nx))
    if source is None:
        source = np.zeros((ny, nx))

    if bc.top is not None:
        T[0, :] = bc.top
    if bc.bottom is not None:
        T[-1, :] = bc.bottom
    if bc.left is not None:
        T[:, 0] = bc.left
    if bc.right is not None:
        T[:, -1] = bc.right

    rx = dy**2 / (2 * (dx**2 + dy**2))
    ry = dx**2 / (2 * (dx**2 + dy**2))
    rs = dx**2 * dy**2 / (2 * k * (dx**2 + dy**2))

    for iteration in range(max_iter):
        max_change = 0.0
        for j in range(1, ny - 1):
            for i in range(1, nx - 1):
                T_old = T[j, i]
                T_left_val = T[j, i - 1]
                T_right_val = T[j, i + 1]
                T_top_val = T[j - 1, i]
                T_bottom_val = T[j + 1, i]

                if bc.left is None and i == 1:
                    T_left_val = T[j, i + 1]
                if bc.right is None and i == nx - 2:
                    T_right_val = T[j, i - 1]
                if bc.top is None and j == 1:
                    T_top_val = T[j + 1, i]
                if bc.bottom is None and j == ny - 2:
                    T_bottom_val = T[j - 1, i]

                T[j, i] = (
                    rx * (T_left_val + T_right_val)
                    + ry * (T_top_val + T_bottom_val)
                    + rs * source[j, i]
                )
                change = abs(T[j, i] - T_old)
                if change > max_change:
                    max_change = change

        if bc.left is None:
            T[:, 0] = T[:, 1]
        if bc.right is None:
            T[:, -1] = T[:, -2]
        if bc.top is None:
            T[0, :] = T[1, :]
        if bc.bottom is None:
            T[-1, :] = T[-2, :]

        if max_change < tol:
            break

    return T


def solve_2d_transient_explicit(T_initial, Lx, Ly, alpha, dt, n_steps, bc, save_every=100):
    """Solve 2D transient conduction with explicit (FTCS) scheme."""
    ny, nx = T_initial.shape
    dx = Lx / (nx - 1)
    dy = Ly / (ny - 1)

    Fo_x = alpha * dt / dx**2
    Fo_y = alpha * dt / dy**2
    if Fo_x + Fo_y > 0.5:
        warnings.warn(
            f"Fourier number Fo_x + Fo_y = {Fo_x + Fo_y:.3f} > 0.5. "
            f"Explicit scheme may be unstable. Reduce dt or increase grid spacing.",
            UserWarning, stacklevel=2,
        )

    T = T_initial.copy()
    history = [T.copy()]

    for step in range(n_steps):
        T_new = T.copy()
        for j in range(1, ny - 1):
            for i in range(1, nx - 1):
                T_new[j, i] = T[j, i] + alpha * dt * (
                    (T[j, i + 1] - 2 * T[j, i] + T[j, i - 1]) / dx**2
                    + (T[j + 1, i] - 2 * T[j, i] + T[j - 1, i]) / dy**2
                )

        if bc.top is not None:
            T_new[0, :] = bc.top
        if bc.bottom is not None:
            T_new[-1, :] = bc.bottom
        if bc.left is not None:
            T_new[:, 0] = bc.left
        if bc.right is not None:
            T_new[:, -1] = bc.right

        T = T_new
        if (step + 1) % save_every == 0 or step == n_steps - 1:
            history.append(T.copy())

    return history
