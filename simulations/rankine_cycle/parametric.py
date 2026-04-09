"""Parametric sweep for Rankine cycle — boiler pressure and superheat temperature."""

from __future__ import annotations
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from simulations.rankine_cycle.rankine import RankineCycle
from thermosim.plotting import apply_style
from thermosim.utils import kelvin_to_celsius, pa_to_mpa


def sweep_boiler_pressure(P_range, T_superheat=773.15, P_condenser=10e3, eta_pump=0.85, eta_turbine=0.87):
    n = len(P_range)
    eta = np.zeros(n)
    w_net = np.zeros(n)
    bwr = np.zeros(n)
    for i, P in enumerate(P_range):
        try:
            result = RankineCycle(fluid="Water", P_boiler=P, T_superheat=T_superheat, P_condenser=P_condenser, eta_pump=eta_pump, eta_turbine=eta_turbine).solve()
            eta[i] = result.thermal_efficiency
            w_net[i] = result.net_work
            bwr[i] = result.back_work_ratio
        except Exception:
            eta[i] = w_net[i] = bwr[i] = np.nan
    return {"P": P_range, "efficiency": eta, "net_work": w_net, "back_work_ratio": bwr}


def sweep_2d(P_range, T_range, P_condenser=10e3, eta_pump=0.85, eta_turbine=0.87):
    nP, nT = len(P_range), len(T_range)
    eta = np.zeros((nP, nT))
    for i, P in enumerate(P_range):
        for j, T in enumerate(T_range):
            try:
                result = RankineCycle(fluid="Water", P_boiler=P, T_superheat=T, P_condenser=P_condenser, eta_pump=eta_pump, eta_turbine=eta_turbine).solve()
                eta[i, j] = result.thermal_efficiency
            except Exception:
                eta[i, j] = np.nan
    return {"P": P_range, "T": T_range, "efficiency": eta}


def main():
    output_dir = Path("simulations/rankine_cycle/outputs")
    output_dir.mkdir(parents=True, exist_ok=True)

    P_range = np.linspace(1e6, 15e6, 30)
    results_1d = sweep_boiler_pressure(P_range, T_superheat=773.15)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    ax1.plot(P_range / 1e6, results_1d["efficiency"] * 100, "o-", color="#e74c3c", linewidth=2)
    apply_style(ax1, title="Thermal Efficiency vs. Boiler Pressure", xlabel="Boiler Pressure [MPa]", ylabel="Thermal Efficiency [%]")
    ax2.plot(P_range / 1e6, results_1d["net_work"] / 1e3, "s-", color="#2980b9", linewidth=2)
    apply_style(ax2, title="Net Work Output vs. Boiler Pressure", xlabel="Boiler Pressure [MPa]", ylabel="Net Work [kJ/kg]")
    fig.tight_layout()
    fig.savefig(output_dir / "parametric_1d.png", dpi=150, bbox_inches="tight")
    print(f"Saved: {output_dir / 'parametric_1d.png'}")
    plt.close(fig)

    P_range_2d = np.linspace(2e6, 14e6, 20)
    T_range_2d = np.linspace(573.15, 873.15, 20)
    results_2d = sweep_2d(P_range_2d, T_range_2d)

    fig, ax = plt.subplots(figsize=(10, 7))
    T_grid, P_grid = np.meshgrid(kelvin_to_celsius(T_range_2d), P_range_2d / 1e6)
    contour = ax.contourf(T_grid, P_grid, results_2d["efficiency"] * 100, levels=20, cmap="RdYlGn")
    fig.colorbar(contour, ax=ax, label="Thermal Efficiency [%]")
    apply_style(ax, title="Rankine Cycle Efficiency Map", xlabel="Superheat Temperature [C]", ylabel="Boiler Pressure [MPa]")
    fig.tight_layout()
    fig.savefig(output_dir / "parametric_2d.png", dpi=150, bbox_inches="tight")
    print(f"Saved: {output_dir / 'parametric_2d.png'}")
    plt.close(fig)


if __name__ == "__main__":
    main()
