"""Entry point for: python -m simulations.rankine_cycle"""

from simulations.rankine_cycle.rankine import RankineCycle
from thermosim.utils import j_to_kj, kelvin_to_celsius, pa_to_mpa


def main():
    print("=" * 60)
    print("  Rankine Cycle Simulation — Single Design Point")
    print("=" * 60)

    cycle = RankineCycle(
        fluid="Water", P_boiler=6e6, T_superheat=773.15,
        P_condenser=10e3, eta_pump=0.85, eta_turbine=0.87,
    )
    result = cycle.solve()

    print("\nOperating Conditions:")
    print(f"  Boiler pressure:    {pa_to_mpa(cycle.P_boiler):.1f} MPa")
    print(f"  Superheat temp:     {kelvin_to_celsius(cycle.T_superheat):.1f} C")
    print(f"  Condenser pressure: {cycle.P_condenser / 1e3:.1f} kPa")
    print(f"  Pump efficiency:    {cycle.eta_pump * 100:.1f}%")
    print(f"  Turbine efficiency: {cycle.eta_turbine * 100:.1f}%")

    print("\nResults:")
    print(f"  Thermal efficiency: {result.thermal_efficiency * 100:.2f}%")
    print(f"  Net work output:    {j_to_kj(result.net_work):.1f} kJ/kg")
    print(f"  Heat input:         {j_to_kj(result.heat_input):.1f} kJ/kg")
    print(f"  Back-work ratio:    {result.back_work_ratio * 100:.2f}%")

    print("\nState Points:")
    labels = ["1 (Pump out)", "2 (Boiler out)", "3 (Turbine out)", "4 (Condenser out)"]
    for label, st in zip(labels, result.states):
        print(
            f"  {label}: T={kelvin_to_celsius(st.T):.1f}C, "
            f"P={pa_to_mpa(st.P):.3f} MPa, "
            f"h={j_to_kj(st.h):.1f} kJ/kg, "
            f"s={st.s / 1e3:.4f} kJ/(kg·K)"
        )

    print("\nEntropy Generation:")
    for comp, ds in result.entropy_generation.items():
        print(f"  {comp}: {ds:.4f} J/(kg·K)")


if __name__ == "__main__":
    main()
