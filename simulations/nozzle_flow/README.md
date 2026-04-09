# Compressible Nozzle Flow Solver

Quasi-1D isentropic flow through a converging-diverging nozzle using area-Mach number relations.

## Usage

```bash
python -m simulations.nozzle_flow
python -m simulations.nozzle_flow --P0 1000000 --T0 400 --A-throat 0.05 --A-exit 0.5
```

## Key Features

- Isentropic area-Mach number relation solver (Brent's method)
- Configurable nozzle geometry with smooth sinusoidal profile
- CLI interface via argparse for parameter specification
- Mach, pressure, and temperature profile visualization
