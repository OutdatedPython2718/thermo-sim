"""Validation of Rankine cycle model against textbook reference data."""

from __future__ import annotations

from dataclasses import dataclass

from simulations.rankine_cycle.rankine import RankineCycle
from thermosim.utils import j_to_kj


@dataclass
class ValidationCase:
    name: str
    P_boiler: float
    T_superheat: float
    P_condenser: float
    eta_pump: float
    eta_turbine: float
    expected_efficiency: float
    expected_net_work: float
    tolerance_efficiency: float
    tolerance_work: float
    source: str


VALIDATION_CASES = [
    ValidationCase(
        name="Ideal Rankine — 6 MPa / 350C / 10 kPa",
        P_boiler=6e6, T_superheat=623.15, P_condenser=10e3,
        eta_pump=1.0, eta_turbine=1.0,
        expected_efficiency=0.363, expected_net_work=1032.0,
        tolerance_efficiency=0.05, tolerance_work=0.05,
        source="Cengel & Boles, Example 10-1 (adapted, CoolProp reference values)",
    ),
    ValidationCase(
        name="Ideal Rankine — 10 MPa / 500C / 10 kPa",
        P_boiler=10e6, T_superheat=773.15, P_condenser=10e3,
        eta_pump=1.0, eta_turbine=1.0,
        expected_efficiency=0.402, expected_net_work=1275.0,
        tolerance_efficiency=0.05, tolerance_work=0.05,
        source="Cengel & Boles, Table 10 (adapted, CoolProp reference values)",
    ),
]


def run_validation():
    results = []
    for case in VALIDATION_CASES:
        result = RankineCycle(
            fluid="Water", P_boiler=case.P_boiler, T_superheat=case.T_superheat,
            P_condenser=case.P_condenser, eta_pump=case.eta_pump, eta_turbine=case.eta_turbine,
        ).solve()
        eta_err = (
            abs(result.thermal_efficiency - case.expected_efficiency) / case.expected_efficiency
        )
        w_err = abs(j_to_kj(result.net_work) - case.expected_net_work) / case.expected_net_work
        results.append({
            "name": case.name, "source": case.source,
            "computed_efficiency": result.thermal_efficiency,
            "expected_efficiency": case.expected_efficiency,
            "efficiency_error_pct": eta_err * 100,
            "efficiency_pass": eta_err < case.tolerance_efficiency,
            "computed_work_kj": j_to_kj(result.net_work),
            "expected_work_kj": case.expected_net_work,
            "work_error_pct": w_err * 100,
            "work_pass": w_err < case.tolerance_work,
        })
    return results


def main():
    print("=" * 70)
    print("  Rankine Cycle Validation — Textbook Comparison")
    print("=" * 70)
    results = run_validation()
    all_pass = True
    for r in results:
        status = "PASS" if (r["efficiency_pass"] and r["work_pass"]) else "FAIL"
        if status == "FAIL":
            all_pass = False
        print(f"\n[{status}] {r['name']}")
        print(f"  Source: {r['source']}")
        print(
            f"  Efficiency: computed={r['computed_efficiency']:.4f}, "
            f"expected={r['expected_efficiency']:.4f}, "
            f"error={r['efficiency_error_pct']:.2f}%"
        )
        print(
            f"  Net work:   computed={r['computed_work_kj']:.1f} kJ/kg, "
            f"expected={r['expected_work_kj']:.1f} kJ/kg, "
            f"error={r['work_error_pct']:.2f}%"
        )
    print(f"\n{'=' * 70}")
    print(f"  Overall: {'ALL PASSED' if all_pass else 'SOME FAILED'}")
    print(f"{'=' * 70}")


if __name__ == "__main__":
    main()
