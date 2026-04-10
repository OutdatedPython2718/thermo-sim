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

## Validation

The model is validated against reference cases adapted from Cengel & Boles *Thermodynamics: An Engineering Approach*:

| Case | Conditions | Expected η | Computed η | Error |
|------|-----------|-----------|-----------|-------|
| Ideal Rankine | 6 MPa / 350°C / 10 kPa | 36.3% | 36.3% | <0.1% |
| Ideal Rankine | 10 MPa / 500°C / 10 kPa | 40.2% | 40.2% | <0.1% |

Both efficiency and net work output are compared, with a 5% relative error tolerance. CoolProp's real-fluid property lookups replace ideal-gas assumptions, so computed values reflect actual steam table data rather than simplified models. Run `python -m simulations.rankine_cycle.validate` to reproduce.

Additional physics checks enforced in the test suite:
- Irreversible cycles (η_pump < 1, η_turbine < 1) always produce lower efficiency than ideal cycles at the same conditions
- Entropy generation is non-negative in every component (2nd law compliance)
- State point pressures match the specified boiler and condenser pressures

## Key Features

- Real-fluid properties via CoolProp (not ideal-gas assumptions)
- Iterative state-point resolution at specified conditions
- 1D and 2D parametric sweeps with efficiency contour maps
- T-S and P-H diagrams with saturation dome overlay
