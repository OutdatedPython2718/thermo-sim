# System-Level Energy Optimization

Optimizes a Rankine power cycle with feedwater preheating using CoolProp real-fluid properties and scipy.optimize.

## Model

The system couples the existing `RankineCycle` model with a feedwater preheating loop:

- **Primary cycle**: Rankine cycle (pump, boiler, turbine, condenser) with CoolProp state-point calculations
- **Heat recovery**: A fraction of the condenser heat rejection is recovered to preheat feedwater before the boiler, reducing heat input
- **Decision variables**: Turbine inlet temperature (550-900 K) and heat recovery fraction (0.0-0.8)
- **Objective**: Minimize specific fuel consumption (heat input per unit net work)
- **Constraint**: Turbine exhaust quality > 0.85 (wet steam erosion limit)

All thermodynamic properties come from CoolProp's IAPWS-IF97 implementation.

## Validation

At zero heat recovery, the model exactly reproduces the standalone `RankineCycle` results (thermal efficiency, net work, and heat input match to machine precision). This self-consistency check runs in CI as a pytest test.

## Scope and Limitations

This model demonstrates constrained optimization of a thermodynamic system using real fluid properties. It is not a production-grade regenerative Rankine cycle simulator. The key differences:

### What this model does

- Recovers a fraction of condenser heat rejection to preheat feedwater before the boiler
- All state points computed via CoolProp's IAPWS-IF97 (real steam properties, not correlations)
- Energy balance is exact: reduced boiler heat input with unchanged turbine/pump work
- Optimization finds the operating point that minimizes fuel consumption subject to a real engineering constraint (turbine exhaust quality > 0.85 to prevent wet steam blade erosion)

### What a production regenerative Rankine model requires

A real power plant feedwater heating system differs from this model in several fundamental ways:

**Bleed steam extraction, not condenser heat recovery.** In practice, steam is extracted from the turbine at one or more intermediate pressures and routed to feedwater heaters. This reduces turbine work output (less mass flow through the low-pressure stages) but improves cycle efficiency by raising the average temperature of heat addition. Our model recovers condenser heat without reducing turbine work — a thermodynamic simplification that overstates the benefit of heat recovery.

**Open and closed feedwater heaters.** Real plants use direct-contact heaters (deaerators, where bleed steam mixes with feedwater) and shell-and-tube closed heaters (where bleed steam condenses on one side). Each heater operates at a different pressure level and requires its own energy and mass balance. A typical power plant has 5-8 feedwater heaters. Our model has a single abstract recovery fraction.

**Multi-stage turbine with extraction points.** Production models split the turbine into high-pressure, intermediate-pressure, and low-pressure sections, each with separate isentropic efficiencies. At each extraction point, a fraction of the steam flow is diverted, and the downstream stages see reduced mass flow. The Baumann rule or similar corrections account for efficiency degradation in wet steam regions.

**Heat exchanger performance parameters.** Real feedwater heaters are characterized by terminal temperature difference (TTD, typically 3-5°C) and drain cooler approach (DCA), which define how close the heater operates to ideal performance. Pressure drops through heaters, pipes, and valves also affect the cycle. Our model assumes ideal heat transfer with no pressure losses.

**Iterative mass flow solution.** With multiple extraction points, the system requires iterative solution for consistent mass flow splits at each bleed point. This is a system of coupled nonlinear equations, typically solved by successive substitution or Newton's method on the full flowsheet.

### Why this model is still valid for its purpose

The optimization methodology (scipy SLSQP with nonlinear constraints), the CoolProp integration, and the feedwater preheating energy balance are all correct. The model demonstrates how to formulate and solve a constrained thermodynamic optimization problem with real fluid properties. The simplification is in the cycle architecture, not in the thermodynamics or the optimization approach.

## Usage

```bash
jupyter notebook simulations/system_optim/system_optimization.ipynb
```

## Key Features

- CoolProp-based thermodynamic model (no simplified correlations)
- scipy.optimize.minimize with SLSQP and real engineering constraints
- Plotly interactive contour plots of the SFC landscape
- Pandas sensitivity analysis tables
