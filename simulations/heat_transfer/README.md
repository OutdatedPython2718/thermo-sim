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
The steady-state solver reproduces the exact linear conduction profile T(x) = T_left + (T_right - T_left) * x / L to machine precision (~10⁻¹² K error).

### Transient Benchmarks

**Semi-infinite solid** (Incropera & DeWitt, Ch. 5, Eq. 5.60): A bar initially at T=0 with one end suddenly set to T=100. The analytical solution T(x,t) = 100 * erfc(x / (2*sqrt(alpha*t))) is compared against the numerical FTCS result at 5 interior positions. Max error < 2%.

**Fourier series finite bar** (Carslaw & Jaeger, Section 3.3): A bar initially at T=0 with T(0)=100, T(L)=0. The exact solution is a Fourier sine series truncated at 50 terms. The numerical solution is compared at moderate time (< 2% error) and long time (< 0.5% error vs. steady-state linear profile).

### Grid Convergence
The 2D Gauss-Seidel solver is run at three resolutions (21x21, 41x41, 81x81). The observed convergence order is computed as p = log(e1/e2) / log(dx1/dx2) and asserted to be between 1.8 and 2.2, confirming 2nd-order spatial accuracy.

All validation tests run in CI via pytest.

## Key Features

- 1D and 2D finite-difference conduction solvers
- Mixed boundary conditions (Dirichlet, Neumann, convective)
- Explicit transient scheme with automatic stability checking
- Heatmaps, cross-section profiles, and convergence plots
