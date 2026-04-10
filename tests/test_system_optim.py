"""Tests for system optimization model."""

import pytest

from simulations.rankine_cycle.rankine import RankineCycle
from simulations.system_optim.model import system_model_with_recovery


class TestSystemModelSelfConsistency:
    """At zero recovery, must match standalone RankineCycle exactly."""

    def test_zero_recovery_matches_standalone(self):
        T_superheat = 773.15
        P_boiler = 6e6
        P_condenser = 10e3

        standalone = RankineCycle(
            fluid="Water",
            P_boiler=P_boiler,
            T_superheat=T_superheat,
            P_condenser=P_condenser,
            eta_pump=0.85,
            eta_turbine=0.87,
        ).solve()

        system_result = system_model_with_recovery(
            T_superheat=T_superheat,
            f_recovery=0.0,
            P_boiler=P_boiler,
            P_condenser=P_condenser,
            eta_pump=0.85,
            eta_turbine=0.87,
        )

        assert system_result["thermal_efficiency"] == pytest.approx(
            standalone.thermal_efficiency, rel=1e-10,
        )
        assert system_result["net_work"] == pytest.approx(
            standalone.net_work, rel=1e-10,
        )
        assert system_result["heat_input"] == pytest.approx(
            standalone.heat_input, rel=1e-10,
        )

    def test_recovery_improves_efficiency(self):
        base = system_model_with_recovery(
            T_superheat=773.15, f_recovery=0.0,
        )
        with_recovery = system_model_with_recovery(
            T_superheat=773.15, f_recovery=0.3,
        )
        assert with_recovery["thermal_efficiency"] > base["thermal_efficiency"]

    def test_exhaust_quality_computed(self):
        result = system_model_with_recovery(
            T_superheat=773.15, f_recovery=0.0,
        )
        assert "exhaust_quality" in result

    def test_all_properties_positive(self):
        result = system_model_with_recovery(
            T_superheat=773.15, f_recovery=0.3,
        )
        assert result["heat_input"] > 0
        assert result["net_work"] > 0
        assert result["thermal_efficiency"] > 0
        assert result["thermal_efficiency"] < 1
