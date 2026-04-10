"""CoolProp-based system model with feedwater preheating.

Models a Rankine cycle where a fraction of the turbine exhaust heat
is recovered to preheat the boiler feedwater, reducing heat input.
"""

from __future__ import annotations

from simulations.rankine_cycle.rankine import RankineCycle


def system_model_with_recovery(
    T_superheat: float = 773.15,
    f_recovery: float = 0.0,
    P_boiler: float = 6e6,
    P_condenser: float = 10e3,
    eta_pump: float = 0.85,
    eta_turbine: float = 0.87,
) -> dict:
    """Compute cycle performance with feedwater preheating.

    Parameters
    ----------
    T_superheat : float
        Turbine inlet temperature [K].
    f_recovery : float
        Fraction of exhaust-to-condenser heat recovered
        for feedwater preheating [0-1].
    P_boiler, P_condenser : float
        Pressures [Pa].
    eta_pump, eta_turbine : float
        Component isentropic efficiencies.

    Returns
    -------
    dict with: thermal_efficiency, net_work, heat_input,
               exhaust_quality, sfc, q_recovered
    """
    cycle = RankineCycle(
        fluid="Water",
        P_boiler=P_boiler,
        T_superheat=T_superheat,
        P_condenser=P_condenser,
        eta_pump=eta_pump,
        eta_turbine=eta_turbine,
    )
    result = cycle.solve()

    state_pump_out = result.states[0]
    state_boiler_out = result.states[1]
    state_turbine_out = result.states[2]
    state_cond_out = result.states[3]

    # Heat available in condenser
    q_rejected = state_turbine_out.h - state_cond_out.h

    # Feedwater preheating
    q_recovered = f_recovery * q_rejected
    h_feedwater_preheated = state_pump_out.h + q_recovered

    # New heat input
    q_in_new = state_boiler_out.h - h_feedwater_preheated

    # Net work unchanged
    w_net = result.net_work

    # New efficiency
    eta_new = w_net / q_in_new if q_in_new > 0 else 0.0

    # SFC in kJ/kWh
    sfc = (q_in_new / max(w_net, 1e-6)) * 3600

    return {
        "thermal_efficiency": eta_new,
        "net_work": w_net,
        "heat_input": q_in_new,
        "exhaust_quality": state_turbine_out.quality,
        "sfc": sfc,
        "q_recovered": q_recovered,
    }
