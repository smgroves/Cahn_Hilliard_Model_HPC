# CPC_geometric_array/

Large parameter sweep over Conserved Phase-field Contact (CPC) geometry simulations. CPC structures are ring-like or crosshair-like domains that model structured biological condensates (e.g., chromosomal cohesin domains). This folder contains the run scripts and job option files for sweeping over CPC ring width, cohesin width, domain size, ε, and α.

---

## Directory Structure

```
CPC_geometric_array/
├── CPC_geometric_run.jl              # Single-job run script (no array; loops internally)
├── CPC_geometric_run.sh              # SLURM submission script for CPC_geometric_run.jl
├── CPC_geometric_run_array.jl        # Array run script; reads CPC_width, cohesin_width, nx from ARGS
├── CPC_geometric_array_*.txt         # Parameter option files for different studies (see below)
└── (no Reports/ here — outputs go to HPC scratch)
```

---

## Run Scripts

### `CPC_geometric_run.jl`
Loops over hard-coded arrays of `CPC_width ∈ {15, 20, 25, 30}` and `cohesin_width ∈ {2, 4, 6, 8}` (grid points). For each combination:
1. Calls `initialize_round_CPC(nx, nx, CPC_width=..., cohesin_width=...)` on a 256² grid
2. Runs `main()` from `../CH_multigrid_solver.jl` with fixed parameters:
   - `dt = 0.1 * (1/256)²`, `gam` set via M=8, `max_it` for `total_time = 0.25`, `tol = 1e-5`
3. Writes output to `/scratch/xpz5km/Cahn_Hilliard_Model/julia_out/critical_radius/CPC_geometry`
4. Appends to `Job_specs.csv` on the HPC home directory

### `CPC_geometric_run_array.jl`
Reads `CPC_width`, `cohesin_width`, `nx` from command-line arguments (`ARGS[1]`, `ARGS[2]`, `ARGS[3]`). Fixed parameters:
- `gam = 0.039117`, `max_it = 20000`, `tol = 1e-5`
- Output: `/scratch/xpz5km/Cahn_Hilliard_Model/julia_out/critical_radius/CPC_geometry_eps_<gam>/`

This script is used when submitting a SLURM array job reading from an options file.

### `CPC_geometric_run.sh`
Single SLURM job (not an array). Requests 16 tasks, 50 GB RAM, 10 hours. Runs `CPC_geometric_run.jl` directly.

---

## Option Files

Each `.txt` file is a parameter list for a batch SLURM array job. Each line specifies arguments for one simulation. The naming convention encodes the study parameters:

| File | Study |
|---|---|
| `CPC_geometric_array_eps_0.0075.txt` | Baseline sweep, ε ≈ 0.0075 |
| `CPC_geometric_array_with_alpha_change_domain_options.txt` | Alpha ≠ 0 with variable domain width |
| `CPC_geometric_array_with_alpha_change_domain_options_eps_0.0075.txt` | Alpha + variable domain, ε ≈ 0.0075 |
| `CPC_geometric_array_with_alpha_change_domain_options_eps_0.0067.txt` | Alpha + variable domain, ε ≈ 0.0067 |
| `CPC_geometric_array_with_alpha_change_domain_options_eps_0.0067_large_widths_t_0.4.txt` | Large CPC widths, run to t = 0.4 |
| `CPC_geometric_array_with_alpha_change_domain_options_eps_0.0075_alpha_-0.5.txt` | α = -0.5 sweep |
| `CPC_geometric_array_eps_0.0075_monse_quant.txt` | Monolayer quantification parameter set |
| `CPC_geometric_array_noisy_cohesin_full_set.txt` | Noisy cohesin initial conditions |
| `CPC_geometric_array_change_domain_options_eps_0.0075_crosshair_localization.txt` | Crosshair geometry + domain change |
| `*_v2.txt / *_v3.txt / *_rerun.txt` | Successive reruns for missing or failed jobs |
| `*_missing_sims.txt / *_redo_more_mem.txt` | Targeted resubmission files |

### Parameter Format (typical)
Lines in the "change domain options" files pass: `CPC_width cohesin_width nx domain_width alpha`

Lines in the simpler eps-only files pass: `CPC_width cohesin_width nx`

---

## Geometry Initialization

`initialize_round_CPC(nx, ny; CPC_width, cohesin_width)` (defined in `CH_multigrid_solver.jl`):
- Creates a ring-shaped CPC domain: inner radius determined by `cohesin_width`, outer radius determined by `CPC_width`
- φ = +1 inside the cohesin region, φ = -1 outside, with a tanh interface
- Models the cohesin ring at the base of a loop domain

`initialize_geometric_CPC(nx, ny; CPC_width, cohesin_width)` (alternate, commented out in most scripts):
- Creates a rectangular/crosshair CPC geometry instead of a circular one

---

## Output

Simulation phi fields are written to HPC scratch storage:
```
/scratch/xpz5km/Cahn_Hilliard_Model/julia_out/critical_radius/CPC_geometry/
```
or for the array variant:
```
/scratch/xpz5km/Cahn_Hilliard_Model/julia_out/critical_radius/CPC_geometry_eps_<gam>/
```

Output filename pattern:
```
phi_<nx>_<max_it>_<tol>__CPC_<CPC_width>_cohesin_<cohesin_width>_eps_<gam>.txt
```

Timing and parameters are appended to `Job_specs.csv` on the HPC login node.
