"""Combined visualization runner for heat transfer simulation."""

from __future__ import annotations

from simulations.heat_transfer.steady_state import main as steady_main
from simulations.heat_transfer.transient import main as transient_main
from simulations.heat_transfer.validate import main as validate_main


def main():
    print("\n" + "=" * 60)
    print("  Generating All Heat Transfer Visualizations")
    print("=" * 60)
    steady_main()
    print()
    transient_main()
    print()
    validate_main()
    print("\n  All heat transfer visualizations complete.")


if __name__ == "__main__":
    main()
