"""Steady-state 2D heat conduction runner with visualization."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from simulations.heat_transfer.conduction import BoundaryCondition, create_grid, solve_2d_steady
from thermosim.plotting import apply_style


def main():
    output_dir = Path("simulations/heat_transfer/outputs")
    output_dir.mkdir(parents=True, exist_ok=True)

    nx, ny = 51, 51
    Lx, Ly = 0.1, 0.1
    k = 50.0

    bc = BoundaryCondition(top=500.0, bottom=300.0, left=300.0, right=400.0)

    print("=" * 60)
    print("  2D Steady-State Heat Conduction")
    print("=" * 60)
    print(f"  Grid: {nx} x {ny}")
    print(f"  Domain: {Lx*100:.0f} cm x {Ly*100:.0f} cm")
    print(f"  k = {k} W/(m·K)")
    print(f"  BCs: top={bc.top}K, bottom={bc.bottom}K, left={bc.left}K, right={bc.right}K")

    grid = create_grid(Lx, Ly, nx, ny)
    T = solve_2d_steady(nx, ny, Lx, Ly, bc, k)

    print(f"\n  T_min = {T.min():.1f} K, T_max = {T.max():.1f} K")
    print(f"  T_center = {T[ny//2, nx//2]:.1f} K")

    fig, ax = plt.subplots(figsize=(8, 7))
    X, Y = np.meshgrid(grid["x"] * 100, grid["y"] * 100)
    contour = ax.contourf(X, Y, T, levels=30, cmap="hot")
    fig.colorbar(contour, ax=ax, label="Temperature [K]")
    contour_lines = ax.contour(X, Y, T, levels=10, colors="white", linewidths=0.5, alpha=0.5)
    ax.clabel(contour_lines, inline=True, fontsize=8, fmt="%.0f")
    apply_style(ax, title="Steady-State Temperature Distribution", xlabel="x [cm]", ylabel="y [cm]")
    fig.tight_layout()
    fig.savefig(output_dir / "steady_state_heatmap.png", dpi=150, bbox_inches="tight")
    print(f"\n  Saved: {output_dir / 'steady_state_heatmap.png'}")
    plt.close(fig)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    mid_j = ny // 2
    ax1.plot(grid["x"] * 100, T[mid_j, :], "o-", color="#e74c3c", markersize=3)
    apply_style(
        ax1,
        title=f"Temperature at y = {grid['y'][mid_j]*100:.1f} cm",
        xlabel="x [cm]",
        ylabel="Temperature [K]",
    )
    mid_i = nx // 2
    ax2.plot(grid["y"] * 100, T[:, mid_i], "s-", color="#2980b9", markersize=3)
    apply_style(
        ax2,
        title=f"Temperature at x = {grid['x'][mid_i]*100:.1f} cm",
        xlabel="y [cm]",
        ylabel="Temperature [K]",
    )
    fig.tight_layout()
    fig.savefig(output_dir / "steady_state_profiles.png", dpi=150, bbox_inches="tight")
    print(f"  Saved: {output_dir / 'steady_state_profiles.png'}")
    plt.close(fig)


if __name__ == "__main__":
    main()
