"""Utility functions for unit conversion, convergence checking, and data export."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np


def celsius_to_kelvin(T_C: float) -> float:
    return T_C + 273.15

def kelvin_to_celsius(T_K: float) -> float:
    return T_K - 273.15

def kpa_to_pa(P_kPa: float) -> float:
    return P_kPa * 1e3

def mpa_to_pa(P_MPa: float) -> float:
    return P_MPa * 1e6

def pa_to_kpa(P_Pa: float) -> float:
    return P_Pa / 1e3

def pa_to_mpa(P_Pa: float) -> float:
    return P_Pa / 1e6

def bar_to_pa(P_bar: float) -> float:
    return P_bar * 1e5

def kj_to_j(E_kJ: float) -> float:
    return E_kJ * 1e3

def j_to_kj(E_J: float) -> float:
    return E_J / 1e3


def check_convergence(values: list[float], tol: float = 1e-6, window: int = 3) -> bool:
    """Check if the last `window` values have converged within tolerance."""
    if len(values) < window:
        return False
    recent = values[-window:]
    return max(recent) - min(recent) < tol


def export_csv(filepath: str | Path, headers: list[str], rows: list[list]) -> None:
    """Export tabular data to CSV."""
    filepath = Path(filepath)
    filepath.parent.mkdir(parents=True, exist_ok=True)
    with open(filepath, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(headers)
        writer.writerows(rows)


def export_json(filepath: str | Path, data: dict) -> None:
    """Export data dict to JSON, handling numpy types."""
    filepath = Path(filepath)
    filepath.parent.mkdir(parents=True, exist_ok=True)

    def default(obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        raise TypeError(f"Not JSON serializable: {type(obj)}")

    with open(filepath, "w") as f:
        json.dump(data, f, indent=2, default=default)
