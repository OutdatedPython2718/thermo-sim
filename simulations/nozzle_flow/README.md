# Compressible Nozzle Flow Solver

Quasi-1D isentropic flow through a converging-diverging nozzle using area-Mach number relations.

## Usage

```bash
python -m simulations.nozzle_flow
python -m simulations.nozzle_flow --P0 1000000 --T0 400 --A-throat 0.05 --A-exit 0.5
```

## Validation

The solver is verified against the analytical isentropic flow relations:

- **Throat condition**: The area-Mach relation requires M = 1.0 at the minimum area (throat). The solver reproduces this within ±0.05 Mach number.
- **Subsonic/supersonic branches**: For a given area ratio A/A* > 1, Brent's method finds both the subsonic (M < 1) and supersonic (M > 1) roots. The test suite verifies the correct branch is selected on each side of the throat.
- **Pressure monotonicity**: In the supersonic diverging section, pressure must decrease monotonically. The test suite confirms P_exit < P_throat for a fully expanded nozzle.
- **Isentropic relations**: Temperature and pressure at every grid point are computed from T/T₀ = (1 + (γ-1)/2 · M²)⁻¹ and P/P₀ = (T/T₀)^(γ/(γ-1)), which are exact for isentropic flow.

## Key Features

- Isentropic area-Mach number relation solver (Brent's method)
- Configurable nozzle geometry with smooth sinusoidal profile
- CLI interface via argparse for parameter specification
- Mach, pressure, and temperature profile visualization
