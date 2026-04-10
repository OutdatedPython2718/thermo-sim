# System-Level Energy Optimization

Jupyter notebook demonstrating scipy.optimize for finding optimal operating points in a combined thermal system (power cycle + waste heat recovery).

## Usage

```bash
jupyter notebook simulations/system_optim/system_optimization.ipynb
```

## Optimization Details

The optimizer minimizes specific fuel consumption (SFC in kJ/kWh) over two decision variables:

- **Turbine inlet temperature** (600–1200 K): Higher temperatures improve Carnot-like efficiency but are bounded by material limits
- **Heat recovery fraction** (0.0–0.8): Preheating the feed with exhaust heat reduces fuel input, but excessive recovery drops exhaust temperature below the 400 K materials constraint

The SLSQP solver enforces the exhaust temperature constraint as a nonlinear inequality (T_exhaust - 400 ≥ 0). The sensitivity analysis table shows how SFC and system efficiency vary across a 5×4 grid of operating points, making the trade-off between temperature and recovery visible.

## Key Features

- Simplified combined cycle thermodynamic model
- scipy.optimize.minimize with SLSQP and nonlinear constraints
- Plotly interactive contour plots of the objective landscape
- Pandas sensitivity analysis tables
