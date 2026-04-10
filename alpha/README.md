# alpha/

Parameter sweep studying how the free-energy asymmetry coefficient **α** affects droplet critical radius behavior. When α ≠ 0, the double-well free energy becomes asymmetric — the two bulk equilibrium phases φ₁ and φ₂ are no longer ±1, and the critical radius for droplet growth or dissolution shifts accordingly.

---

## Directory Structure

```
alpha/
├── CH_multigrid_solver_with_alpha_v2.jl     # Local copy of the alpha-enabled multigrid solver
├── julia_run_critical_radius_alpha_v3.jl    # Run script: droplet critical radius with alpha
├── run_julia_multigrid_alpha_v3.sh          # SLURM batch script (array job)
├── options_alpha_+0.5.txt                   # Parameter list for α = +0.5 runs
├── options_alpha_+0.5_v2.txt                # Updated parameter list for α = +0.5
├── options_alpha_-0.5.txt                   # Parameter list for α = -0.5 runs
├── options_alpha_-0.5_v2.txt / v3 / v4      # Successive revisions of α = -0.5 parameters
├── level_set_plot_alpha_v3.m                # MATLAB: plot droplet radius vs. time (alpha runs)
├── level_set_radius_array_alpha_v3.m        # MATLAB: radius computation over parameter arrays
├── level_set_radius_plot_array_alpha_v3.sh  # Shell: MATLAB batch launcher for radius arrays
├── level_set_plot_array/                    # Subdirectory with per-parameter-set level set plots
├── droplet_with_atanh_interface/            # Simulations using atanh-profile droplet IC
├── flat_interface/                          # Simulations with flat interface IC
└── Reports/                                 # Simulation output organized by Slurm job ID
```

---

## Solver

`CH_multigrid_solver_with_alpha_v2.jl` is a local copy of `../CH_multigrid_solver_with_alpha_change_domain.jl`. The modified free energy derivative is:

```
f'(φ) = φ³ - φ + α
```

For non-zero α, the two equilibrium phases φ₁ and φ₂ (roots of f'(φ) = 0) shift symmetrically. The initialization function `initialization_from_function_v2()` accepts these phase values directly to construct the correct tanh droplet profile.

---

## Run Script: `julia_run_critical_radius_alpha_v3.jl`

Accepts command-line arguments: `R0 gam total_time alpha nx`

For each job:
1. Selects precomputed equilibrium parameters (β, φ₁, φ₂) for the given α:
   - α = -0.5: β = 0.6794, φ₁ = -1.2075, φ₂ = 0.8742
   - α = -0.1: β = 0.7024, φ₁ = -1.0350, φ₂ = 0.9683
   - α = +0.1: β = 0.7059, φ₁ = -0.9683, φ₂ = 1.0350
   - α = +0.5: β = 0.6794, φ₁ = -0.8742, φ₂ = 1.2075
2. Initializes a circular droplet with radius R0 using the asymmetric tanh profile
3. Runs the multigrid solver via `main_w_alpha()`
4. Writes output to `/project/g_bme-janeslab/SarahG/julia_out/critical_radius_alpha_updated/`
5. Appends timing and parameters to `/home/xpz5km/Cahn_Hilliard_Model/Job_specs.csv`

Output filename pattern:
```
phi_<nx>_<max_it>_<tol>__R0_<R0>_eps_<gam>_alpha_<alpha>.txt
```

---

## SLURM Job: `run_julia_multigrid_alpha_v3.sh`

| Setting | Value |
|---|---|
| Nodes | 1 |
| Tasks per node | 16 |
| Memory | 200 GB |
| Wall time | 10 hours |
| Array size | 1–78 (reads from `options_alpha_+0.5_v2.txt`) |
| SLURM output | `./Reports/%A/output.%J.out` |

Array task ID maps to a line in the options file. Each line specifies: `R0 gam total_time alpha nx`.

---

## Option Files

Each line represents one simulation job: `R0  gam  total_time  alpha  nx`

Example from `options_alpha_-0.5.txt`:
```
0.004 0.0018761 10 -0.5 256
0.005 0.0018761 10 -0.5 256
...
```

The `gam` values correspond to different values of M (number of interface mesh points) at given `nx`:
- `gam ≈ 0.0019` → M = 4 at nx = 256 (thinner interface)
- `gam ≈ 0.0075` → M = 8 at nx = 256 (standard interface)
- `gam ≈ 0.015`  → M = 4 at nx = 128

---

## Post-Processing

### `level_set_radius_array_alpha_v3.m`
Loops over all (R0, gam, alpha) combinations, reads the phi output, and computes the droplet radius at each saved timepoint using the zero level set (φ = 0 contour). Results are aggregated for plotting.

### `level_set_plot_alpha_v3.m`
Plots radius vs. time for each run, overlaying different alpha values to visualize how asymmetry shifts the critical radius threshold (the radius below which droplets dissolve vs. grow).

### `level_set_radius_plot_array_alpha_v3.sh`
Calls MATLAB in batch mode to run the radius analysis scripts on the cluster.

---

## Subdirectories

### `droplet_with_atanh_interface/`
Simulations initialized with a droplet using an atanh interface profile (the proper equilibrium profile for the double-well potential), rather than a sharp step. Ensures the IC is consistent with the steady-state interface shape.

### `flat_interface/`
Simulations with a flat interface IC, used to study how the interface relaxes to its equilibrium profile when α ≠ 0.
