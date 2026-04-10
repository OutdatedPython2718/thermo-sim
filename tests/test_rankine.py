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
