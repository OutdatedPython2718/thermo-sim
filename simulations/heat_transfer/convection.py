"""Convective boundary condition handling and heat transfer coefficient correlations."""

from __future__ import annotations


def htc_flat_plate_forced(Re, Pr, k_fluid, L):
    """Average HTC for forced convection over a flat plate (laminar). Re < 5e5."""
    Nu = 0.664 * Re**0.5 * Pr**(1 / 3)
    return Nu * k_fluid / L


def htc_natural_convection_vertical(Ra, Pr, k_fluid, L):
    """Average HTC for natural convection on a vertical plate (Churchill-Chu)."""
    Nu = (0.825 + 0.387 * Ra**(1 / 6) / (1 + (0.492 / Pr)**(9 / 16))**(8 / 27))**2
    return Nu * k_fluid / L


def apply_convective_bc(T, boundary, h_conv, T_inf, k_solid, dx):
    """Apply convective (Robin) BC to one edge of a 2D temperature field.
    Uses ghost-node approach: -k dT/dn = h(T_surface - T_inf)
    """
    Bi = h_conv * dx / k_solid
    if boundary == "top":
        T[0, :] = (T[1, :] + Bi * T_inf) / (1 + Bi)
    elif boundary == "bottom":
        T[-1, :] = (T[-2, :] + Bi * T_inf) / (1 + Bi)
    elif boundary == "left":
        T[:, 0] = (T[:, 1] + Bi * T_inf) / (1 + Bi)
    elif boundary == "right":
        T[:, -1] = (T[:, -2] + Bi * T_inf) / (1 + Bi)
    else:
        raise ValueError(f"Unknown boundary: {boundary}")
    return T
