"""Transient 2D heat conduction runner."""

from __future__ import annotations
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from simulations.heat_transfer.conduction import BoundaryCondition, create_grid, solve_2d_transient_explicit
from thermosim.plotting import apply_style


def main():
    output_dir = Path("simulations/heat_transfer/outputs")
    output_dir.mkdir(parents=True, exist_ok=True)

    nx, ny = 31, 31
    Lx, Ly = 0.1, 0.1
    alpha = 1.17e-5
    dx = Lx / (nx - 1)
    dt_max = 0.25 * dx**2 / alpha
    dt = 0.8 * dt_max
    n_steps = 3000
    save_every = 300

    bc = BoundaryCondition(top=500.0, bottom=300.0, left=300.0, right=300.0)
    T_initial = np.full((ny, nx), 300.0)
    T_initial[0, :] = 500.0

    print("=" * 60)
    print("  2D Transient Heat Conduction (Explicit FTCS)")
    print("=" * 60)
    print(f"  Grid: {nx} x {ny}")
    print(f"  dt = {dt:.6f} s, n_steps = {n_steps}")
    print(f"  Fo_x = {alpha * dt / dx**2:.4f}")
    print(f"  Total time = {dt * n_steps:.3f} s")

    history = solve_2d_transient_explicit(
        T_initial=T_initial, Lx=Lx, Ly=Ly, alpha=alpha,
        dt=dt, n_steps=n_steps, bc=bc, save_every=save_every,
    )

    n_snapshots = min(len(history), 6)
    indices = np.linspace(0, len(history) - 1, n_snapshots, dtype=int)
    grid = create_grid(Lx, Ly, nx, ny)

    fig, axes = plt.subplots(2, 3, figsize=(16, 10))
    axes = axes.flatten()
    X, Y = np.meshgrid(grid["x"] * 100, grid["y"] * 100)
    vmin = T_initial.min()
    vmax = bc.top

    for idx, (ax, snap_i) in enumerate(zip(axes[:n_snapshots], indices)):
        T_snap = history[snap_i]
        time_s = snap_i * save_every * dt if snap_i > 0 else 0
        contour = ax.contourf(X, Y, T_snap, levels=20, cmap="hot", vmin=vmin, vmax=vmax)
        apply_style(ax, title=f"t = {time_s:.3f} s", xlabel="x [cm]", ylabel="y [cm]")

    fig.colorbar(contour, ax=axes.tolist(), label="Temperature [K]", shrink=0.8)
    fig.suptitle("Transient Temperature Evolution", fontsize=16, fontweight="bold", y=1.01)
    fig.tight_layout()
    fig.savefig(output_dir / "transient_snapshots.png", dpi=150, bbox_inches="tight")
    print(f"\n  Saved: {output_dir / 'transient_snapshots.png'}")
    plt.close(fig)

    center_temps = [h[ny // 2, nx // 2] for h in history]
    times = [i * save_every * dt for i in range(len(history))]
    times[0] = 0.0

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(times, center_temps, "o-", color="#e74c3c", markersize=4)
    apply_style(ax, title="Center Temperature vs. Time", xlabel="Time [s]", ylabel="Temperature [K]")
    fig.tight_layout()
    fig.savefig(output_dir / "transient_center_temp.png", dpi=150, bbox_inches="tight")
    print(f"  Saved: {output_dir / 'transient_center_temp.png'}")
    plt.close(fig)


if __name__ == "__main__":
    main()
