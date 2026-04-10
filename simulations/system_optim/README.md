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

## Usage

```bash
jupyter notebook simulations/system_optim/system_optimization.ipynb
```

## Key Features

- CoolProp-based thermodynamic model (no simplified correlations)
- scipy.optimize.minimize with SLSQP and real engineering constraints
- Plotly interactive contour plots of the SFC landscape
- Pandas sensitivity analysis tables
