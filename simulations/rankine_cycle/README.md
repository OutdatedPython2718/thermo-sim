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

State points and cycle efficiency are validated against NIST/IAPWS-IF97 steam table reference values (source: [NIST Chemistry Webbook](https://webbook.nist.gov/chemistry/fluid/)).

| State Point | Condition | NIST h [kJ/kg] | NIST s [kJ/(kg·K)] |
|-------------|-----------|----------------|---------------------|
| Sat. liquid | 10 kPa | 191.83 | 0.6493 |
| Superheated | 6 MPa, 350°C | 3043.0 | 6.3335 |
| Superheated | 10 MPa, 500°C | 3373.7 | 6.5966 |

Each state point's enthalpy and entropy computed by CoolProp is compared against the NIST value within 0.5% relative error. Overall cycle efficiency is verified against hand calculations using the NIST enthalpy values:

| Case | Conditions | NIST η | Computed η | Tolerance |
|------|-----------|--------|-----------|-----------|
| Ideal Rankine | 6 MPa / 350°C / 10 kPa | 36.27% | ~36.3% | < 2% |
| Ideal Rankine | 10 MPa / 500°C / 10 kPa | 40.18% | ~40.2% | < 2% |

These tests run in CI via pytest. Run `python -m simulations.rankine_cycle.validate` to reproduce.

## Key Features

- Real-fluid properties via CoolProp (not ideal-gas assumptions)
- Iterative state-point resolution at specified conditions
- 1D and 2D parametric sweeps with efficiency contour maps
- T-S and P-H diagrams with saturation dome overlay
