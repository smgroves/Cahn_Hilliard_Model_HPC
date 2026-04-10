# Cahn-Hilliard Model HPC

High-performance computing implementations of solvers for the **Cahn-Hilliard equation**, a fourth-order nonlinear PDE used to model phase separation and coarsening dynamics in binary mixtures. This repository is used to study phenomena such as spinodal decomposition, droplet critical radius evolution, and Conserved Phase-field Contact (CPC) geometry—with relevance to biological structures such as chromosomal condensates.

**Author:** Sarah Groves

**Note: All code for execution was written by Sarah Groves. Claude was used to generate READMEs for general use.**
---

## Mathematical Background

The Cahn-Hilliard equation governs the evolution of an order parameter φ ∈ [-1, 1] (e.g., relative concentration):

```
∂φ/∂t = ∇²μ
μ = f'(φ) - ε²∇²φ
```

where:
- μ is the chemical potential
- f'(φ) = φ³ - φ is the derivative of the double-well free energy (optionally modified with an asymmetry parameter α)
- ε is the interfacial thickness parameter, set via `gam = M * h / (2√2 * atanh(0.9))` where `M` is the number of grid points spanning the interface and `h = 1/nx` is the grid spacing

---

## Repository Structure

```
.
├── CH_multigrid_solver.jl                    # Core multigrid solver (Neumann BCs)
├── CH_multigrid_solver_periodic_bc.jl        # Core multigrid solver (periodic BCs)
├── CH_multigrid_solver_with_alpha_change_domain.jl  # Solver with alpha asymmetry + domain options
├── CHplotting_function.m                     # MATLAB: generate AVI movies from phi output
├── CHplotting_function_minutes.m             # Variant with time in minutes
├── level_set_plot.m                          # MATLAB: droplet radius via level-set (φ=0)
├── level_set_plot_alt_IC.m                   # Variant for alternate initial conditions
├── level_set_radius.m                        # Compute level-set radius at each timeframe
├── level_set_radius_array.m                  # Radius computation over parameter arrays
├── level_set_radius_array_alt_IC.m           # Radius array computation with alt IC
├── level_set_radius_multiple_droplets.m      # Radius tracking for multi-droplet runs
├── plot_IC_heatmap.m                         # Visualize initial condition heatmap
├── plot_single_timepoint_heatmap.m           # Visualize a single timepoint heatmap
├── make_plots.sh                             # Shell: batch plotting launcher
├── make_plots_critical_radius.sh             # Shell: plots for critical radius study
├── make_plots_crosshair.sh                   # Shell: crosshair cross-section plots
├── make_plots_Figure_6_supp.sh               # Shell: manuscript supplementary figure
├── level_set_radius_plot.sh                  # Shell: run level_set_radius.m via MATLAB CLI
├── level_set_radius_plot_array.sh            # Shell: run level_set_radius_array.m
├── level_set_radius_plot_array_alt_IC.sh     # Shell: alt IC array run
├── Job_specs.csv                             # Master log of simulation parameters and timing
├── frame.png                                 # Sample phi field visualization
├── spinodal_decomp/                          # Reference implementations and cross-language comparison
├── evaluating_phi/                           # Studies on phi behavior and alpha effects
├── alpha/                                    # Parameter sweep over alpha (free energy asymmetry)
├── CPC_geometric_array/                      # CPC geometry array simulations
├── CPC_1024_geometry/                        # High-resolution (1024²) CPC geometry run
├── run_julia_multigrid_critical_radius/      # Critical radius study with multiple configurations
└── Reports/                                  # Simulation output organized by Slurm job ID
```

---

## Core Solvers

All three top-level `.jl` files follow the same architecture. Each is self-contained and `include`d by run scripts in the subfolders.

### `CH_multigrid_solver.jl` — Neumann BCs
The primary solver. Key functions:

| Function | Description |
|---|---|
| `laplace(a, nx, ny)` | 5-point finite difference Laplacian with homogeneous Neumann (zero-flux) boundary conditions |
| `source(c_old, nx, ny, dt)` | Assembles right-hand side for the discretized system |
| `relax(c, mu, su, sw, ...)` | Damped Jacobi relaxation on the coupled (φ, μ) 2×2 block system |
| `restrict_ch(a, nx, ny)` | Restriction operator: injection from fine to coarse grid |
| `prolong_ch(a, nx, ny)` | Prolongation operator: bilinear interpolation from coarse to fine |
| `vcycle(c, mu, su, sw, ...)` | Recursive V-cycle multigrid algorithm |
| `error2(c, mu, su, sw, ...)` | Computes residual norm for convergence check |
| `initialization(...)` | Creates initial conditions: `"random"`, `"droplet"` (tanh profile), `"spinodal"`, geometric CPC |
| `initialize_round_CPC(...)` | Circular CPC domain initialization |
| `initialize_geometric_CPC(...)` | Rectangular/crosshair CPC domain initialization |
| `cahn(phi, nx, tol, outdir, ...)` | Main time-stepping loop; calls `vcycle` each step until residual < `tol` |
| `main(phi, nx, tol, outdir, ...)` | Entry point; handles I/O, output directories, timing |

Output files are written to `outdir` and named by parameter suffix (e.g., `phi_128_500_1.0e-6__R0_0.1_eps_0.015.txt`). Each file contains the full phi field concatenated over all saved timesteps.

### `CH_multigrid_solver_periodic_bc.jl` — Periodic BCs
Identical structure to the Neumann solver; `laplace()` uses periodic index wraparound instead of zero-flux conditions. Used for spinodal decomposition and unbounded-domain studies.

### `CH_multigrid_solver_with_alpha_change_domain.jl` — Alpha + Domain Options
Extends the Neumann solver with:
- A free-energy asymmetry parameter `alpha` that shifts f'(φ) = φ³ - φ + α
- Support for non-unit domain widths (used in CPC geometry studies where physical domain proportions vary)

---

## MATLAB Post-Processing

| Script | Purpose |
|---|---|
| `CHplotting_function.m` | Reads concatenated phi text files, reshapes into frames, generates an AVI movie with a red-blue colormap (range [-1, 1]) |
| `CHplotting_function_minutes.m` | Same but labels time axis in minutes |
| `level_set_radius.m` | Computes droplet radius at each timeframe by finding the zero level set of φ |
| `level_set_radius_array.m` | Loops over parameter combinations (R0, eps) and collects radii |
| `level_set_radius_multiple_droplets.m` | Handles multi-droplet runs where multiple zero contours exist |
| `level_set_plot.m` | Plots radius vs. time for a single run |
| `plot_IC_heatmap.m` | Displays the initial φ field as a heatmap |
| `plot_single_timepoint_heatmap.m` | Extracts and displays a single saved timepoint |

---

## Job Tracking

`Job_specs.csv` is the master log appended by all run scripts. Columns include: simulation type, language, grid size (`nx`), `gam` (ε), `dt`, tolerance, max iterations, CH iterations, and wall-clock time.

---

## HPC Environment

Simulations are run on a SLURM cluster. Key paths written by run scripts:
- **Scratch output:** `/scratch/xpz5km/Cahn_Hilliard_Model/julia_out/`
- **Project storage:** `/project/g_bme-janeslab/SarahG/julia_out/`
- **Job specs log:** `/home/xpz5km/Cahn_Hilliard_Model/Job_specs.csv`

Julia module used: `julia/1.9.2`

---

## Subfolders

| Folder | Summary |
|---|---|
| [`spinodal_decomp/`](spinodal_decomp/README.md) | Multi-language solver implementations and cross-language accuracy comparison |
| [`evaluating_phi/`](evaluating_phi/README.md) | Analysis of φ behavior under varying alpha and geometry |
| [`alpha/`](alpha/README.md) | Parameter sweep over the free-energy asymmetry parameter α |
| [`CPC_geometric_array/`](CPC_geometric_array/README.md) | Large array of CPC geometry simulations across domain/width parameters |
| [`CPC_1024_geometry/`](CPC_1024_geometry/README.md) | High-resolution 1024² CPC geometry run for manuscript |
| [`run_julia_multigrid_critical_radius/`](run_julia_multigrid_critical_radius/README.md) | Critical radius study with droplet ICs and multiple grid/BC configurations |
| [`Reports/`](Reports/) | Raw simulation outputs organized by Slurm job ID (not documented here) |
