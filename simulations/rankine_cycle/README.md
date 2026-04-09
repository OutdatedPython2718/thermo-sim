# Rankine Cycle Simulation

A superheated Rankine cycle model using CoolProp real-fluid properties for water/steam. Computes thermal efficiency, net work output, back-work ratio, and entropy generation with configurable component efficiencies.

## Physics

The model simulates four processes:
1. **Pump** (1→2): Adiabatic compression of saturated liquid to boiler pressure
2. **Boiler** (2→3): Constant-pressure heat addition to superheated vapor
3. **Turbine** (3→4): Adiabatic expansion to condenser pressure
4. **Condenser** (4→1): Constant-pressure heat rejection to saturated liquid

Isentropic efficiencies for pump and turbine account for irreversibilities.

## Usage

```bash
# Single design point
python -m simulations.rankine_cycle

# Parametric sweep (boiler pressure and superheat temperature)
python -m simulations.rankine_cycle.parametric

# Validation against textbook data
python -m simulations.rankine_cycle.validate

# Generate T-S, P-H diagrams and entropy generation chart
python -m simulations.rankine_cycle.visualize
```

## Key Features

- Real-fluid properties via CoolProp (not ideal-gas assumptions)
- Iterative state-point resolution at specified conditions
- 1D and 2D parametric sweeps with efficiency contour maps
- Validation against Cengel & Boles textbook reference values
- Publication-quality T-S and P-H diagrams with saturation dome
