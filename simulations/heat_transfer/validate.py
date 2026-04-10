"""Validation of heat transfer solvers against analytical solutions."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from simulations.heat_transfer.conduction import BoundaryCondition, solve_1d_steady, solve_2d_steady
from thermosim.plotting import apply_style


def analytical_1d_conduction(x, T_left, T_right, L):
    return T_left + (T_right - T_left) * x / L


def validate_1d():
    L = 1.0
    T_left, T_right = 100.0, 500.0
    nx = 51
    T_numerical = solve_1d_steady(nx=nx, L=L, T_left=T_left, T_right=T_right)
    x = np.linspace(0, L, nx)
    T_analytical = analytical_1d_conduction(x, T_left, T_right, L)
    error = np.max(np.abs(T_numerical - T_analytical))
    return {
        "x": x,
        "T_numerical": T_numerical,
        "T_analytical": T_analytical,
        "max_error": error,
        "pass": error < 0.01,
    }


def grid_convergence_study():
    bc = BoundaryCondition(top=400.0, bottom=200.0, left=300.0, right=100.0)
    resolutions = [11, 21, 41, 81, 161]
    center_temps = []
    dx_values = []
    for n in resolutions:
        T = solve_2d_steady(nx=n, ny=n, Lx=1.0, Ly=1.0, bc=bc, k=1.0, tol=1e-8, max_iter=50000)
        center_temps.append(T[n // 2, n // 2])
        dx_values.append(1.0 / (n - 1))
    T_exact_est = center_temps[-1] + (center_temps[-1] - center_temps[-2]) / 3
    errors = [abs(tc - T_exact_est) for tc in center_temps]
    return {
        "resolutions": resolutions,
        "dx_values": dx_values,
        "center_temps": center_temps,
        "errors": errors,
        "T_exact_est": T_exact_est,
    }


def main():
    output_dir = Path("simulations/heat_transfer/outputs")
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("  Heat Transfer Validation & Grid Convergence")
    print("=" * 60)

    print("\n--- 1D Steady Conduction Validation ---")
    result_1d = validate_1d()
    status = "PASS" if result_1d["pass"] else "FAIL"
    print(f"  [{status}] Max error vs. analytical: {result_1d['max_error']:.2e} K")

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(
        result_1d["x"], result_1d["T_analytical"],
        "--", color="#2c3e50", linewidth=2, label="Analytical",
    )
    ax.plot(
        result_1d["x"], result_1d["T_numerical"],
        "o", color="#e74c3c", markersize=4, label="Numerical",
    )
    apply_style(
        ax,
        title="1D Conduction: Numerical vs. Analytical",
        xlabel="x [m]",
        ylabel="Temperature [K]",
    )
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_dir / "validation_1d.png", dpi=150, bbox_inches="tight")
    print(f"  Saved: {output_dir / 'validation_1d.png'}")
    plt.close(fig)

    print("\n--- 2D Grid Convergence Study ---")
    gc = grid_convergence_study()
    print(f"  Richardson-extrapolated center T = {gc['T_exact_est']:.4f} K")
    for n, tc, err in zip(gc["resolutions"], gc["center_temps"], gc["errors"]):
        print(f"  {n:>4d}x{n:<4d}: T_center = {tc:.4f} K, error = {err:.2e}")

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.loglog(gc["dx_values"], gc["errors"], "o-", color="#e74c3c", linewidth=2, markersize=8)
    dx_ref = np.array(gc["dx_values"])
    scale = gc["errors"][0] / dx_ref[0]**2
    ax.loglog(dx_ref, scale * dx_ref**2, "--", color="#7f8c8d", label="2nd order ref.")
    apply_style(
        ax,
        title="Grid Convergence Study",
        xlabel="Grid Spacing dx [m]",
        ylabel="Error vs. Extrapolated [K]",
    )
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_dir / "grid_convergence.png", dpi=150, bbox_inches="tight")
    print(f"  Saved: {output_dir / 'grid_convergence.png'}")
    plt.close(fig)


if __name__ == "__main__":
    main()
