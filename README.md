# Thermo-Sim

**Thermodynamic simulation toolkit** — cycle analysis, heat transfer modeling, and system optimization in Python.

![Python](https://img.shields.io/badge/Python-3.11+-blue?logo=python&logoColor=white)
![NumPy](https://img.shields.io/badge/NumPy-1.24+-013243?logo=numpy)
![SciPy](https://img.shields.io/badge/SciPy-1.10+-8CAAE6?logo=scipy)
![CoolProp](https://img.shields.io/badge/CoolProp-6.6+-orange)
![matplotlib](https://img.shields.io/badge/matplotlib-3.7+-11557c)
![License](https://img.shields.io/badge/License-MIT-green)

## Overview

A collection of thermodynamic simulations demonstrating real-fluid property modeling, numerical methods, and engineering analysis. Built as a modular Python package with shared utilities, validated solvers, and clear visualizations.

```
thermo-sim/
├── thermosim/           # Shared library: fluid properties, solvers, plotting
├── simulations/
│   ├── rankine_cycle/   # Superheated Rankine cycle with CoolProp
│   ├── heat_transfer/   # 2D finite-difference conduction (steady + transient)
│   ├── nozzle_flow/     # Quasi-1D compressible nozzle flow
│   └── system_optim/    # Combined cycle optimization (scipy.optimize)
├── notebooks/           # Exploratory analysis with Plotly
└── tests/               # pytest unit and integration tests
```

## Gallery

| Rankine Cycle T-S Diagram | 2D Temperature Field | Nozzle Mach Profile |
|:-:|:-:|:-:|
| ![T-S](docs/gallery/ts_diagram.png) | ![Heatmap](docs/gallery/steady_state_heatmap.png) | ![Nozzle](docs/gallery/nozzle_profiles.png) |

## Quick Start

```bash
# Clone and install
git clone https://github.com/USERNAME/thermo-sim.git
cd thermo-sim
pip install -e ".[dev]"

# Run a simulation
python -m simulations.rankine_cycle
python -m simulations.heat_transfer.steady_state
python -m simulations.nozzle_flow
```

## Simulations

### [Rankine Cycle](simulations/rankine_cycle/)
Superheated Rankine cycle with CoolProp real-fluid properties, parametric sweeps over boiler pressure and superheat temperature, and state-point validation against NIST/IAPWS-IF97 steam tables (enthalpy and entropy within 0.5%, cycle efficiency within 2% of hand calculation from NIST reference values).

### [Heat Transfer](simulations/heat_transfer/)
1D/2D finite-difference conduction with steady-state (Gauss-Seidel) and transient (explicit FTCS) solvers. Validated against three analytical benchmarks: exact 1D linear profile (machine precision), semi-infinite solid erfc solution (< 2% error), and Fourier series finite bar (< 2% error). Grid convergence study confirms 2nd-order spatial accuracy with computed convergence order p ≈ 2.0.

### [Nozzle Flow](simulations/nozzle_flow/)
Quasi-1D compressible flow through a converging-diverging nozzle. Solves the isentropic area-Mach relation using Brent's root-finding method, then computes pressure, temperature, and density from isentropic relations. Verified by checking Mach = 1.0 at the throat and confirming pressure/temperature ratios match the analytical isentropic relations at the exit plane.

### [System Optimization](simulations/system_optim/)
Rankine cycle with CoolProp-based feedwater preheating, optimized via scipy.optimize with a turbine exhaust quality constraint. Self-consistency verified: zero-recovery case exactly matches standalone Rankine cycle results.

### [Exploratory Analysis](notebooks/)
Working fluid comparison (Water, R-134a, CO2) with Pandas and Plotly — saturation pressure curves, latent heat profiles, critical point properties, and overlaid T-s saturation envelopes.

## Testing

```bash
pytest tests/ -v
```

## Docker

```bash
docker build -t thermosim .
docker run thermosim
docker run thermosim python -m simulations.heat_transfer.steady_state
```

## Tech Stack

| Category | Tools |
|----------|-------|
| Scientific Computing | NumPy, SciPy, CoolProp |
| Visualization | matplotlib, Plotly |
| Data Analysis | Pandas, Jupyter |
| Infrastructure | Git, Docker, GitHub Actions CI/CD |
| Code Quality | pytest, ruff |

## License

MIT
