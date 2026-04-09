"""Rankine cycle simulation with CoolProp real-fluid properties.

Models a simple superheated Rankine cycle: pump -> boiler -> turbine -> condenser.
"""

from __future__ import annotations
from dataclasses import dataclass
from thermosim.fluids import FluidState, fluid_state


@dataclass
class CycleResult:
    """Results from a Rankine cycle calculation."""
    states: list[FluidState]  # [pump_out, boiler_out, turbine_out, condenser_out]
    thermal_efficiency: float
    net_work: float          # J/kg
    heat_input: float        # J/kg
    heat_rejected: float     # J/kg
    back_work_ratio: float
    entropy_generation: dict[str, float]


class RankineCycle:
    def __init__(self, fluid="Water", P_boiler=6e6, T_superheat=773.15, P_condenser=10e3, eta_pump=0.85, eta_turbine=0.87):
        self.fluid = fluid
        self.P_boiler = P_boiler
        self.T_superheat = T_superheat
        self.P_condenser = P_condenser
        self.eta_pump = eta_pump
        self.eta_turbine = eta_turbine

    def solve(self) -> CycleResult:
        # State 3 (condenser outlet): saturated liquid at condenser pressure
        state_cond_out = fluid_state(self.fluid, P=self.P_condenser, Q=0.0)

        # Pump: isentropic compression
        state_pump_ideal = fluid_state(self.fluid, P=self.P_boiler, s=state_cond_out.s)
        w_pump_ideal = state_pump_ideal.h - state_cond_out.h
        w_pump_actual = w_pump_ideal / self.eta_pump
        h_pump_out = state_cond_out.h + w_pump_actual
        state_pump_out = fluid_state(self.fluid, P=self.P_boiler, h=h_pump_out)

        # Boiler outlet: superheated vapor
        state_boiler_out = fluid_state(self.fluid, P=self.P_boiler, T=self.T_superheat)

        # Turbine: isentropic expansion
        state_turb_ideal = fluid_state(self.fluid, P=self.P_condenser, s=state_boiler_out.s)
        w_turb_ideal = state_boiler_out.h - state_turb_ideal.h
        w_turb_actual = w_turb_ideal * self.eta_turbine
        h_turb_out = state_boiler_out.h - w_turb_actual
        state_turb_out = fluid_state(self.fluid, P=self.P_condenser, h=h_turb_out)

        # Energy balance
        q_in = state_boiler_out.h - state_pump_out.h
        q_out = state_turb_out.h - state_cond_out.h
        w_net = w_turb_actual - w_pump_actual
        eta_thermal = w_net / q_in
        bwr = w_pump_actual / w_turb_actual

        entropy_gen = {
            "pump": state_pump_out.s - state_cond_out.s,
            "turbine": state_turb_out.s - state_boiler_out.s,
        }

        return CycleResult(
            states=[state_pump_out, state_boiler_out, state_turb_out, state_cond_out],
            thermal_efficiency=eta_thermal, net_work=w_net, heat_input=q_in,
            heat_rejected=q_out, back_work_ratio=bwr, entropy_generation=entropy_gen,
        )
