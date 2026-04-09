"""Quasi-1D compressible nozzle flow solver."""

from __future__ import annotations
import argparse
from dataclasses import dataclass
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq
from thermosim.plotting import apply_style


@dataclass
class NozzleGeometry:
    L: float
    A_inlet: float
    A_throat: float
    A_exit: float

    def area(self, x):
        x_norm = x / self.L
        x_throat = 0.4
        return np.where(
            x_norm <= x_throat,
            self.A_inlet + (self.A_throat - self.A_inlet) * (1 - np.cos(np.pi * x_norm / x_throat)) / 2,
            self.A_throat + (self.A_exit - self.A_throat) * (1 - np.cos(np.pi * (x_norm - x_throat) / (1 - x_throat))) / 2,
        )


def isentropic_mach_from_area_ratio(A_ratio, gamma=1.4, supersonic=False):
    gp1 = gamma + 1
    gm1 = gamma - 1

    def equation(M):
        return (1 / M) * ((2 / gp1) * (1 + gm1 / 2 * M**2)) ** (gp1 / (2 * gm1)) - A_ratio

    if A_ratio < 1.0 + 1e-10:
        return 1.0
    if supersonic:
        return brentq(equation, 1.0, 20.0)
    else:
        return brentq(equation, 0.01, 1.0)


def solve_nozzle_flow(geom, P0=500e3, T0=300.0, gamma=1.4, R=287.0, n_points=200):
    gm1 = gamma - 1
    x = np.linspace(0, geom.L, n_points)
    A = geom.area(x)
    A_star = geom.A_throat
    Mach = np.zeros(n_points)
    throat_idx = np.argmin(A)

    for i in range(n_points):
        A_ratio = A[i] / A_star
        supersonic = i > throat_idx
        try:
            Mach[i] = isentropic_mach_from_area_ratio(A_ratio, gamma, supersonic)
        except ValueError:
            Mach[i] = Mach[i - 1] if i > 0 else 0.01

    T = T0 / (1 + gm1 / 2 * Mach**2)
    P = P0 * (T / T0) ** (gamma / gm1)
    rho = P / (R * T)
    V = Mach * np.sqrt(gamma * R * T)
    return {"x": x, "A": A, "Mach": Mach, "P": P, "T": T, "rho": rho, "V": V}


def main():
    parser = argparse.ArgumentParser(description="Quasi-1D compressible nozzle flow solver")
    parser.add_argument("--P0", type=float, default=500e3, help="Stagnation pressure [Pa]")
    parser.add_argument("--T0", type=float, default=300.0, help="Stagnation temperature [K]")
    parser.add_argument("--gamma", type=float, default=1.4, help="Specific heat ratio")
    parser.add_argument("--A-inlet", type=float, default=0.2, help="Inlet area [m2]")
    parser.add_argument("--A-throat", type=float, default=0.1, help="Throat area [m2]")
    parser.add_argument("--A-exit", type=float, default=0.3, help="Exit area [m2]")
    parser.add_argument("--length", type=float, default=1.0, help="Nozzle length [m]")
    parser.add_argument("-n", type=int, default=200, help="Number of grid points")
    args = parser.parse_args()

    geom = NozzleGeometry(L=args.length, A_inlet=args.A_inlet, A_throat=args.A_throat, A_exit=args.A_exit)
    result = solve_nozzle_flow(geom=geom, P0=args.P0, T0=args.T0, gamma=args.gamma, n_points=args.n)

    output_dir = Path("simulations/nozzle_flow/outputs")
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("  Quasi-1D Compressible Nozzle Flow")
    print("=" * 60)
    print(f"  P0 = {args.P0/1e3:.1f} kPa, T0 = {args.T0:.1f} K, gamma = {args.gamma}")
    print(f"  Exit Mach = {result['Mach'][-1]:.3f}")
    print(f"  Exit P = {result['P'][-1]/1e3:.1f} kPa")
    print(f"  Exit T = {result['T'][-1]:.1f} K")

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 12), sharex=True)
    ax1.fill_between(result["x"], result["A"], alpha=0.3, color="#3498db")
    ax1.plot(result["x"], result["A"], color="#2980b9", linewidth=2)
    apply_style(ax1, title="Nozzle Area Distribution", ylabel="Area [m2]")

    ax2.plot(result["x"], result["Mach"], color="#e74c3c", linewidth=2)
    ax2.axhline(y=1.0, color="#7f8c8d", linestyle="--", alpha=0.5, label="M = 1")
    apply_style(ax2, title="Mach Number Distribution", ylabel="Mach Number")
    ax2.legend()

    ax3.plot(result["x"], result["P"] / args.P0, color="#27ae60", linewidth=2)
    apply_style(ax3, title="Pressure Ratio Distribution", xlabel="Axial Position [m]", ylabel="P/P0")

    fig.tight_layout()
    fig.savefig(output_dir / "nozzle_profiles.png", dpi=150, bbox_inches="tight")
    print(f"\n  Saved: {output_dir / 'nozzle_profiles.png'}")
    plt.close(fig)


if __name__ == "__main__":
    main()
