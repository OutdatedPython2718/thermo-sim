"""Tests for Rankine cycle simulation."""

import pytest

from simulations.rankine_cycle.rankine import CycleResult, RankineCycle


class TestRankineCycle:
    def test_basic_cycle_runs(self):
        cycle = RankineCycle(
            fluid="Water", P_boiler=6e6, T_superheat=600.0 + 273.15,
            P_condenser=10e3, eta_pump=0.85, eta_turbine=0.87,
        )
        result = cycle.solve()
        assert isinstance(result, CycleResult)
        assert result.thermal_efficiency > 0
        assert result.thermal_efficiency < 1
        assert result.net_work > 0
        assert result.heat_input > 0
        assert result.back_work_ratio > 0
        assert result.back_work_ratio < 1
        assert len(result.states) == 4

    def test_efficiency_in_expected_range(self):
        cycle = RankineCycle(
            fluid="Water", P_boiler=6e6, T_superheat=773.15,
            P_condenser=10e3, eta_pump=1.0, eta_turbine=1.0,
        )
        result = cycle.solve()
        assert 0.30 < result.thermal_efficiency < 0.45

    def test_irreversible_less_efficient_than_ideal(self):
        ideal = RankineCycle(
            fluid="Water", P_boiler=6e6, T_superheat=773.15,
            P_condenser=10e3, eta_pump=1.0, eta_turbine=1.0,
        ).solve()
        real = RankineCycle(
            fluid="Water", P_boiler=6e6, T_superheat=773.15,
            P_condenser=10e3, eta_pump=0.85, eta_turbine=0.87,
        ).solve()
        assert real.thermal_efficiency < ideal.thermal_efficiency

    def test_entropy_generation_positive(self):
        cycle = RankineCycle(
            fluid="Water", P_boiler=4e6, T_superheat=673.15,
            P_condenser=10e3, eta_pump=0.85, eta_turbine=0.85,
        )
        result = cycle.solve()
        for name, ds in result.entropy_generation.items():
            assert ds >= -1e-6, f"Negative entropy generation in {name}"

    def test_states_have_correct_pressures(self):
        P_boil = 8e6
        P_cond = 15e3
        cycle = RankineCycle(
            fluid="Water", P_boiler=P_boil, T_superheat=773.15,
            P_condenser=P_cond, eta_pump=0.9, eta_turbine=0.9,
        )
        result = cycle.solve()
        assert result.states[0].P == pytest.approx(P_boil, rel=0.01)
        assert result.states[1].P == pytest.approx(P_boil, rel=0.01)
        assert result.states[2].P == pytest.approx(P_cond, rel=0.01)
        assert result.states[3].P == pytest.approx(P_cond, rel=0.01)


class TestRankineNISTValidation:
    """Validate CoolProp (IAPWS-IF97) output against NIST/IAPWS-IF97 reference values.

    Reference values are sourced independently from the NIST Chemistry Webbook
    (https://webbook.nist.gov/chemistry/fluid/) and standard steam tables.
    Tolerance is 0.5% to accommodate minor differences between IAPWS-IF97
    implementations and NIST Webbook rounding.
    """

    def test_condenser_outlet_matches_nist(self):
        """Saturated liquid at 10 kPa: NIST h_f=191.83 kJ/kg, s_f=0.6493 kJ/(kg·K)."""
        from thermosim.fluids import fluid_state

        nist_h = 191830.0   # J/kg  (191.83 kJ/kg)
        nist_s = 649.3      # J/(kg·K)  (0.6493 kJ/(kg·K))

        state = fluid_state("Water", P=10e3, Q=0.0)

        assert state.h == pytest.approx(nist_h, rel=0.005), (
            f"h at sat liq 10 kPa: CoolProp={state.h:.2f} J/kg, NIST={nist_h} J/kg"
        )
        assert state.s == pytest.approx(nist_s, rel=0.005), (
            f"s at sat liq 10 kPa: CoolProp={state.s:.4f} J/(kg·K), NIST={nist_s} J/(kg·K)"
        )

    def test_case1_boiler_outlet_matches_nist(self):
        """Superheated steam at 6 MPa, 350°C: NIST h=3043.0 kJ/kg, s=6.3335 kJ/(kg·K)."""
        from thermosim.fluids import fluid_state

        nist_h = 3043000.0  # J/kg  (3043.0 kJ/kg)
        nist_s = 6333.5     # J/(kg·K)  (6.3335 kJ/(kg·K))

        state = fluid_state("Water", P=6e6, T=623.15)

        assert state.h == pytest.approx(nist_h, rel=0.005), (
            f"h at 6 MPa/350°C: CoolProp={state.h:.2f} J/kg, NIST={nist_h} J/kg"
        )
        assert state.s == pytest.approx(nist_s, rel=0.005), (
            f"s at 6 MPa/350°C: CoolProp={state.s:.4f} J/(kg·K), NIST={nist_s} J/(kg·K)"
        )

    def test_case2_boiler_outlet_matches_nist(self):
        """Superheated steam at 10 MPa, 500°C: NIST h=3373.7 kJ/kg, s=6.5966 kJ/(kg·K)."""
        from thermosim.fluids import fluid_state

        nist_h = 3373700.0  # J/kg  (3373.7 kJ/kg)
        nist_s = 6596.6     # J/(kg·K)  (6.5966 kJ/(kg·K))

        state = fluid_state("Water", P=10e6, T=773.15)

        assert state.h == pytest.approx(nist_h, rel=0.005), (
            f"h at 10 MPa/500°C: CoolProp={state.h:.2f} J/kg, NIST={nist_h} J/kg"
        )
        assert state.s == pytest.approx(nist_s, rel=0.005), (
            f"s at 10 MPa/500°C: CoolProp={state.s:.4f} J/(kg·K), NIST={nist_s} J/(kg·K)"
        )

    def test_case1_cycle_efficiency_from_nist(self):
        """Ideal Rankine at 6 MPa/350°C/10 kPa: expected η ≈ 0.3627.

        Hand calc from NIST steam tables:
          w_turbine = h1 - h2s = 3043.0 - 2005.15 = 1037.85 kJ/kg  (isentropic expansion)
          w_pump    = h4 - h3  = 197.88 - 191.83  = 6.05 kJ/kg      (pump work)
          q_in      = h1 - h4  = 3043.0 - 197.88  = 2845.12 kJ/kg
          η = (w_turbine - w_pump) / q_in = 1031.8 / 2845.12 = 0.3627
        Source: NIST Webbook (https://webbook.nist.gov/chemistry/fluid/)
        """
        nist_eta = 0.3627
        cycle = RankineCycle(
            fluid="Water", P_boiler=6e6, T_superheat=623.15,
            P_condenser=10e3, eta_pump=1.0, eta_turbine=1.0,
        )
        result = cycle.solve()
        assert result.thermal_efficiency == pytest.approx(nist_eta, rel=0.02), (
            f"Case 1 efficiency: computed={result.thermal_efficiency:.4f}, NIST={nist_eta}"
        )

    def test_case2_cycle_efficiency_from_nist(self):
        """Ideal Rankine at 10 MPa/500°C/10 kPa: expected η ≈ 0.4018.

        Hand calc from NIST steam tables:
          w_turbine = h1 - h2s = 3373.7 - 2089.08 = 1284.62 kJ/kg  (isentropic expansion)
          w_pump    = h4 - h3  = 201.92 - 191.83  = 10.09 kJ/kg     (pump work)
          q_in      = h1 - h4  = 3373.7 - 201.92  = 3171.78 kJ/kg
          η = (w_turbine - w_pump) / q_in = 1274.5 / 3171.78 = 0.4018
        Source: NIST Webbook (https://webbook.nist.gov/chemistry/fluid/)
        """
        nist_eta = 0.4018
        cycle = RankineCycle(
            fluid="Water", P_boiler=10e6, T_superheat=773.15,
            P_condenser=10e3, eta_pump=1.0, eta_turbine=1.0,
        )
        result = cycle.solve()
        assert result.thermal_efficiency == pytest.approx(nist_eta, rel=0.02), (
            f"Case 2 efficiency: computed={result.thermal_efficiency:.4f}, NIST={nist_eta}"
        )
