"""Thermodynamic diagram generators.

Provides T-S and P-H diagram plotting with saturation domes,
process paths, and consistent styling.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from thermosim.fluids import saturation_curve

COLORS = {
    "dome": "#2c3e50",
    "process": "#e74c3c",
    "points": "#2980b9",
    "grid": "#ecf0f1",
    "fill": "#3498db",
}


def apply_style(ax, title="", xlabel="", ylabel=""):
    """Apply consistent styling to an axes."""
    ax.set_title(title, fontsize=14, fontweight="bold", pad=12)
    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.tick_params(labelsize=10)
    ax.grid(True, alpha=0.3, color=COLORS["grid"])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def plot_ts_diagram(fluid, states=None, figsize=(10, 7)):
    """Plot a Temperature-Entropy diagram with saturation dome."""
    fig, ax = plt.subplots(figsize=figsize)

    sat = saturation_curve(fluid, n_points=200)
    s_dome = np.concatenate([sat["s_f"], sat["s_g"][::-1]]) / 1e3
    T_dome = np.concatenate([sat["T"], sat["T"][::-1]])

    ax.plot(sat["s_f"] / 1e3, sat["T"], color=COLORS["dome"], linewidth=2, label="Saturation dome")
    ax.plot(sat["s_g"] / 1e3, sat["T"], color=COLORS["dome"], linewidth=2)
    ax.fill(s_dome, T_dome, alpha=0.08, color=COLORS["fill"])

    if states:
        s_vals = [st.s / 1e3 for st in states]
        T_vals = [st.T for st in states]
        ax.plot(s_vals, T_vals, "o-", color=COLORS["process"], linewidth=2,
                markersize=8, markerfacecolor=COLORS["points"], markeredgecolor="white",
                markeredgewidth=1.5, label="Process path", zorder=5)
        for i, (s, T) in enumerate(zip(s_vals, T_vals)):
            ax.annotate(
                f"  {i + 1}", (s, T), fontsize=10, fontweight="bold", color=COLORS["points"]
            )

    apply_style(
        ax, title=f"T-s Diagram — {fluid}", xlabel="Entropy [kJ/(kg·K)]", ylabel="Temperature [K]"
    )
    ax.legend(fontsize=10)
    fig.tight_layout()
    return fig, ax


def plot_ph_diagram(fluid, states=None, figsize=(10, 7)):
    """Plot a Pressure-Enthalpy diagram with saturation dome."""
    fig, ax = plt.subplots(figsize=figsize)

    sat = saturation_curve(fluid, n_points=200)
    ax.plot(
        sat["h_f"] / 1e3, sat["P"] / 1e6, color=COLORS["dome"], linewidth=2, label="Saturation dome"
    )
    ax.plot(sat["h_g"] / 1e3, sat["P"] / 1e6, color=COLORS["dome"], linewidth=2)

    h_dome = np.concatenate([sat["h_f"], sat["h_g"][::-1]]) / 1e3
    P_dome = np.concatenate([sat["P"], sat["P"][::-1]]) / 1e6
    ax.fill(h_dome, P_dome, alpha=0.08, color=COLORS["fill"])

    if states:
        h_vals = [st.h / 1e3 for st in states]
        P_vals = [st.P / 1e6 for st in states]
        ax.plot(h_vals, P_vals, "o-", color=COLORS["process"], linewidth=2,
                markersize=8, markerfacecolor=COLORS["points"], markeredgecolor="white",
                markeredgewidth=1.5, label="Process path", zorder=5)
        for i, (h, P) in enumerate(zip(h_vals, P_vals)):
            ax.annotate(
                f"  {i + 1}", (h, P), fontsize=10, fontweight="bold", color=COLORS["points"]
            )

    ax.set_yscale("log")
    apply_style(
        ax, title=f"P-h Diagram — {fluid}", xlabel="Enthalpy [kJ/kg]", ylabel="Pressure [MPa]"
    )
    ax.legend(fontsize=10)
    fig.tight_layout()
    return fig, ax
