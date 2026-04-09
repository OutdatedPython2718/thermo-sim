# Heat Transfer & Thermal Management Model

Finite-difference solvers for 1D and 2D heat conduction with steady-state and transient analysis. Includes convective boundary conditions, grid convergence studies, and validation against analytical solutions.

## Physics

Models heat conduction governed by:
- **Steady-state**: nabla(k*nabla(T)) + q = 0
- **Transient**: rho*cp * dT/dt = nabla(k*nabla(T)) + q

Boundary conditions: Dirichlet (fixed temperature), Neumann (insulated), and Robin (convective).

## Numerical Methods

- **Steady-state**: Gauss-Seidel iterative solver and direct tridiagonal solve (1D)
- **Transient**: Forward-Time Central-Space (FTCS) explicit scheme with Fourier number stability check
- **Validation**: Grid convergence study with Richardson extrapolation

## Usage

```bash
# Steady-state 2D plate
python -m simulations.heat_transfer.steady_state

# Transient simulation with snapshots
python -m simulations.heat_transfer.transient

# Validation and grid convergence study
python -m simulations.heat_transfer.validate

# Generate all visualizations
python -m simulations.heat_transfer.visualize
```

## Key Features

- 1D and 2D finite-difference conduction solvers
- Mixed boundary conditions (Dirichlet, Neumann, convective)
- Explicit transient scheme with automatic stability checking
- Grid convergence study demonstrating 2nd-order spatial accuracy
- Validation against analytical solutions with quantified error
- Publication-quality heatmaps, cross-section profiles, and convergence plots
