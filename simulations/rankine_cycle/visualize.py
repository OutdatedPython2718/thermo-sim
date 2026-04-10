"""Visualization dashboard for Rankine cycle results."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt

from simulations.rankine_cycle.rankine import RankineCycle
from thermosim.plotting import COLORS, apply_style, plot_ph_diagram, plot_ts_diagram


def plot_entropy_generation(result, output_path=None):
    fig, ax = plt.subplots(figsize=(8, 5))
    components = list(result.entropy_generation.keys())
    values = [result.entropy_generation[c] for c in components]
    bars = ax.bar(
        components, values,
        color=[COLORS["process"], COLORS["points"]], edgecolor="white", linewidth=1.5,
    )
    for bar, val in zip(bars, values):
        ax.text(
            bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.5,
            f"{val:.2f}", ha="center", fontsize=11, fontweight="bold",
        )
    apply_style(
        ax, title="Entropy Generation by Component",
        xlabel="Component", ylabel="Entropy Generation [J/(kg·K)]",
    )
    fig.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=150, bbox_inches="tight")
        print(f"Saved: {output_path}")
    return fig, ax


def main():
    output_dir = Path("simulations/rankine_cycle/outputs")
    output_dir.mkdir(parents=True, exist_ok=True)

    cycle = RankineCycle(
        fluid="Water", P_boiler=6e6, T_superheat=773.15,
        P_condenser=10e3, eta_pump=0.85, eta_turbine=0.87,
    )
    result = cycle.solve()
    plot_states = result.states + [result.states[0]]

    fig_ts, _ = plot_ts_diagram("Water", states=plot_states)
    fig_ts.savefig(output_dir / "ts_diagram.png", dpi=150, bbox_inches="tight")
    print(f"Saved: {output_dir / 'ts_diagram.png'}")
    plt.close(fig_ts)

    fig_ph, _ = plot_ph_diagram("Water", states=plot_states)
    fig_ph.savefig(output_dir / "ph_diagram.png", dpi=150, bbox_inches="tight")
    print(f"Saved: {output_dir / 'ph_diagram.png'}")
    plt.close(fig_ph)

    plot_entropy_generation(result, output_dir / "entropy_generation.png")
    plt.close("all")
    print("\nAll Rankine cycle visualizations generated.")


if __name__ == "__main__":
    main()
