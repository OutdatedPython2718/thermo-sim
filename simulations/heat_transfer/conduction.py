"""Finite-difference conduction solvers for 1D and 2D heat transfer."""

from __future__ import annotations

import warnings
from dataclasses import dataclass

import numpy as np
from scipy.sparse import diags as sp_diags, kron as sp_kron, eye as sp_eye
from scipy.sparse.linalg import spsolve


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
    """Solve 2D steady-state conduction using a sparse direct solver.

    Assembles the second-order central-difference system for the 2D
    conduction equation ``k ∇²T + Q = 0`` and solves it exactly via sparse
    LU factorisation.  Dirichlet conditions are applied where the
    corresponding ``BoundaryCondition`` field is not ``None``; a zero-flux
    (Neumann) mirror condition is used where the field is ``None``.

    To maximise floating-point accuracy the linear system is solved for
    ``T' = T - T_shift``, where ``T_shift`` is the mean of all non-None
    boundary values.  This centring reduces round-off when the boundary
    values are large relative to their differences.

    The ``tol`` and ``max_iter`` parameters are accepted for API compatibility
    but are not used; the direct solver always returns the exact FD solution
    to machine precision.
    """
    dx = Lx / (nx - 1)
    dy = Ly / (ny - 1)
    if source is None:
        source = np.zeros((ny, nx))

    # Collect the active (Dirichlet) boundary values and compute a shift so
    # that the shifted BCs are centred near zero for better round-off.
    bc_vals = [v for v in (bc.top, bc.bottom, bc.left, bc.right) if v is not None]
    T_shift = float(np.mean(bc_vals)) if bc_vals else 0.0

    # Shifted boundary values (None stays None)
    s_top = (bc.top - T_shift) if bc.top is not None else None
    s_bot = (bc.bottom - T_shift) if bc.bottom is not None else None
    s_lft = (bc.left - T_shift) if bc.left is not None else None
    s_rgt = (bc.right - T_shift) if bc.right is not None else None

    ax = 1.0 / dx**2
    ay = 1.0 / dy**2

    ni_x = nx - 2  # number of interior columns
    ni_y = ny - 2  # number of interior rows
    n_int = ni_x * ni_y

    def flat(j, i):
        """Map interior 1-based indices (j, i) to flat index."""
        return (j - 1) * ni_x + (i - 1)

    # Build sparse system in LIL format for efficient row-wise fill.
    from scipy.sparse import lil_matrix
    A = lil_matrix((n_int, n_int))
    b = np.zeros(n_int)

    for j in range(1, ny - 1):
        for i in range(1, nx - 1):
            row = flat(j, i)
            A[row, row] = -2.0 * (ax + ay)

            # Left neighbour
            if i - 1 == 0:
                if s_lft is not None:
                    b[row] -= ax * s_lft
                else:                        # Neumann: mirror
                    A[row, row] += ax
            else:
                A[row, flat(j, i - 1)] += ax

            # Right neighbour
            if i + 1 == nx - 1:
                if s_rgt is not None:
                    b[row] -= ax * s_rgt
                else:
                    A[row, row] += ax
            else:
                A[row, flat(j, i + 1)] += ax

            # Top neighbour (j-1 in array index)
            if j - 1 == 0:
                if s_top is not None:
                    b[row] -= ay * s_top
                else:
                    A[row, row] += ay
            else:
                A[row, flat(j - 1, i)] += ay

            # Bottom neighbour
            if j + 1 == ny - 1:
                if s_bot is not None:
                    b[row] -= ay * s_bot
                else:
                    A[row, row] += ay
            else:
                A[row, flat(j + 1, i)] += ay

            # Source term: k∇²T = -Q  →  ∇²T' = -Q/k  (shift doesn't affect ∇²)
            b[row] -= source[j, i] / k

    T_prime = spsolve(A.tocsr(), b)

    # Reconstruct full temperature field
    T = np.full((ny, nx), T_shift)
    if bc.top is not None:
        T[0, :] = bc.top
    if bc.bottom is not None:
        T[-1, :] = bc.bottom
    if bc.left is not None:
        T[:, 0] = bc.left
    if bc.right is not None:
        T[:, -1] = bc.right
    for j in range(1, ny - 1):
        for i in range(1, nx - 1):
            T[j, i] = T_prime[flat(j, i)] + T_shift

    # Neumann ghost rows/cols
    if bc.left is None:
        T[:, 0] = T[:, 1]
    if bc.right is None:
        T[:, -1] = T[:, -2]
    if bc.top is None:
        T[0, :] = T[1, :]
    if bc.bottom is None:
        T[-1, :] = T[-2, :]

    return T


def solve_1d_transient_explicit(
    T_initial, L, alpha, dt, n_steps,
    T_left, T_right, save_every=100,
):
    """Solve 1D transient conduction with explicit (FTCS) scheme.

    Parameters
    ----------
    T_initial : ndarray, shape (nx,)
        Initial temperature distribution.
    L : float
        Domain length [m].
    alpha : float
        Thermal diffusivity [m^2/s].
    dt : float
        Time step [s].
    n_steps : int
        Number of time steps.
    T_left, T_right : float
        Dirichlet boundary temperatures [K].
    save_every : int
        Save snapshot every N steps.

    Returns
    -------
    history : list of ndarray
        Temperature snapshots.
    """
    nx = len(T_initial)
    dx = L / (nx - 1)

    Fo = alpha * dt / dx**2
    if Fo > 0.5:
        warnings.warn(
            f"Fourier number Fo = {Fo:.3f} > 0.5. "
            f"Explicit scheme may be unstable.",
            UserWarning, stacklevel=2,
        )

    T = T_initial.copy()
    T[0] = T_left
    T[-1] = T_right
    history = [T.copy()]

    for step in range(n_steps):
        T_new = T.copy()
        for i in range(1, nx - 1):
            T_new[i] = T[i] + Fo * (T[i + 1] - 2 * T[i] + T[i - 1])
        T_new[0] = T_left
        T_new[-1] = T_right
        T = T_new

        if (step + 1) % save_every == 0 or step == n_steps - 1:
            history.append(T.copy())

    return history


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
