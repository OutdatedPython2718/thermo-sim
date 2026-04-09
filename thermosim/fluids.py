"""CoolProp wrappers for thermodynamic fluid property lookups.

Provides a FluidState dataclass and helper functions that replace
ideal-gas assumptions with real-fluid property calculations.
"""

from dataclasses import dataclass

import CoolProp.CoolProp as CP
import numpy as np


@dataclass(frozen=True)
class FluidState:
    """Thermodynamic state point for a fluid."""
    fluid: str
    T: float      # Temperature [K]
    P: float      # Pressure [Pa]
    h: float      # Specific enthalpy [J/kg]
    s: float      # Specific entropy [J/(kg·K)]
    rho: float    # Density [kg/m³]
    quality: float # Vapor quality [-1 if superheated/subcooled, 0-1 if two-phase]


_FLUID_MAP = {
    "Water": "Water",
    "R134a": "R134a",
    "CO2": "CO2",
    "Ammonia": "Ammonia",
    "R410A": "R410A",
}

_INPUT_KEYS = {
    "T": "T", "P": "P", "h": "H", "s": "S", "Q": "Q", "rho": "D",
}


def fluid_state(fluid: str, **kwargs) -> FluidState:
    """Compute a full thermodynamic state from two independent properties.

    Parameters
    ----------
    fluid : str
        Fluid name (e.g., "Water", "R134a", "CO2").
    **kwargs : float
        Exactly two of: T [K], P [Pa], h [J/kg], s [J/(kg·K)], Q [-], rho [kg/m³].
    """
    if fluid not in _FLUID_MAP:
        raise ValueError(f"Unknown fluid '{fluid}'. Supported: {list(_FLUID_MAP.keys())}")
    if len(kwargs) != 2:
        raise ValueError(f"Exactly two independent properties required, got {len(kwargs)}: {list(kwargs.keys())}")

    cp_fluid = _FLUID_MAP[fluid]
    keys = list(kwargs.keys())
    vals = list(kwargs.values())
    k1, k2 = _INPUT_KEYS[keys[0]], _INPUT_KEYS[keys[1]]

    T = CP.PropsSI("T", k1, vals[0], k2, vals[1], cp_fluid)
    P = CP.PropsSI("P", k1, vals[0], k2, vals[1], cp_fluid)
    h = CP.PropsSI("H", k1, vals[0], k2, vals[1], cp_fluid)
    s = CP.PropsSI("S", k1, vals[0], k2, vals[1], cp_fluid)
    rho = CP.PropsSI("D", k1, vals[0], k2, vals[1], cp_fluid)
    quality = CP.PropsSI("Q", k1, vals[0], k2, vals[1], cp_fluid)

    return FluidState(fluid=fluid, T=T, P=P, h=h, s=s, rho=rho, quality=quality)


def saturation_curve(fluid: str, n_points: int = 200) -> dict[str, np.ndarray]:
    """Compute saturation dome data for plotting T-S and P-H diagrams."""
    if fluid not in _FLUID_MAP:
        raise ValueError(f"Unknown fluid '{fluid}'. Supported: {list(_FLUID_MAP.keys())}")

    cp_fluid = _FLUID_MAP[fluid]
    T_min = CP.PropsSI("Tmin", cp_fluid) + 1.0
    T_crit = CP.PropsSI("Tcrit", cp_fluid) - 1.0

    temps = np.linspace(T_min, T_crit, n_points)
    s_f = np.zeros(n_points)
    s_g = np.zeros(n_points)
    h_f = np.zeros(n_points)
    h_g = np.zeros(n_points)
    pressures = np.zeros(n_points)

    for i, T in enumerate(temps):
        s_f[i] = CP.PropsSI("S", "T", T, "Q", 0, cp_fluid)
        s_g[i] = CP.PropsSI("S", "T", T, "Q", 1, cp_fluid)
        h_f[i] = CP.PropsSI("H", "T", T, "Q", 0, cp_fluid)
        h_g[i] = CP.PropsSI("H", "T", T, "Q", 1, cp_fluid)
        pressures[i] = CP.PropsSI("P", "T", T, "Q", 0, cp_fluid)

    return {"T": temps, "P": pressures, "s_f": s_f, "s_g": s_g, "h_f": h_f, "h_g": h_g}
