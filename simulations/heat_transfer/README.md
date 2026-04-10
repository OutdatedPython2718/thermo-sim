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

## Validation

### 1D Analytical Comparison

The 1D steady-state solver is validated against the exact analytical solution for conduction through a bar with fixed-temperature boundaries and no internal heat source: T(x) = T_left + (T_right - T_left) * x / L. The tridiagonal matrix solve reproduces this linear profile to near machine precision (~10⁻¹² K max error), confirming the discretization and solver are implemented correctly.

### 2D Grid Convergence Study

The 2D Gauss-Seidel solver is run at five resolutions (11×11, 21×21, 41×41, 81×81, 161×161) on a plate with four different boundary temperatures. The center-point temperature is tracked at each resolution to measure convergence:

1. **Richardson extrapolation** on the two finest grids estimates the "exact" solution
2. Error vs. grid spacing is plotted on a log-log scale
3. The slope is compared against a 2nd-order reference line

This demonstrates mesh independence — the solution converges at the expected O(Δx²) rate for a second-order central-difference scheme. This is standard Verification & Validation (V&V) practice in computational engineering.

### Transient Stability

The explicit FTCS scheme automatically checks the Fourier number (Fo = α·Δt/Δx²) before time-stepping and raises a warning if Fo_x + Fo_y > 0.5, which would violate the CFL stability condition and produce oscillating, non-physical results.

Run `python -m simulations.heat_transfer.validate` to reproduce all validation results.

## Key Features

- 1D and 2D finite-difference conduction solvers
- Mixed boundary conditions (Dirichlet, Neumann, convective)
- Explicit transient scheme with automatic stability checking
- Heatmaps, cross-section profiles, and convergence plots
