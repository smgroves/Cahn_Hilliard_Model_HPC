# run_julia_multigrid_critical_radius/

Studies the **critical radius** phenomenon in the Cahn-Hilliard model: a droplet smaller than the critical radius R* will dissolve, while a droplet larger than R* will grow. This folder contains multiple run scripts for different configurations (grid sizes, boundary conditions, initial condition types, and multi-droplet systems), along with the corresponding SLURM option files and post-processing scripts.

---

## Directory Structure

```
run_julia_multigrid_critical_radius/
├── julia_run_critical_radius.jl                  # Standard single-droplet run (nx=128)
├── julia_run_critical_radius_256.jl              # Single-droplet run at nx=256
├── julia_run_critical_radius_alt_IC.jl           # Alternate initial condition
├── julia_run_critical_radius_offcenter.jl        # Off-center droplet placement
├── julia_run_critical_radius_periodic.jl         # Periodic boundary conditions
├── julia_run_critical_radius_two_drops.jl        # Two-droplet system (nx=512)
├── julia_run_critical_radius_two_drops_256_doubledomain.jl  # Two drops, doubled domain
├── run_julia_multigrid.sh                        # SLURM array script (standard)
├── run_julia_multigrid_256.sh                    # SLURM array script (256 grid)
├── run_julia_multigrid_alt_IC.sh                 # SLURM array script (alt IC)
├── options.txt / options_v2.txt ... v9.txt       # Parameter lists for standard runs
├── options_256.txt / options_256_v2.txt          # Parameter lists for 256-grid runs
├── options_alt_IC_128.txt / options_alt_IC_256.txt  # Parameter lists for alt IC
├── options_large_eps.txt / options_large_eps_v2.txt # Parameter lists for large ε
├── options_periodic_check.txt                    # Parameter list for periodic BC check
├── level_set_radius_multiple_droplets.m          # MATLAB: track radius for multi-droplet runs
├── plot_timepoints.sh                            # Shell: plot selected timepoints
└── plot_timepoints_critical_radius.py            # Python: generate timepoint plots
```

Output is written to `Reports/<JobID>/output.<JobID>.out` by SLURM.

---

## Run Scripts

### `julia_run_critical_radius.jl` — Standard Run (nx=128)

Accepts: `R0  M  total_time` as command-line arguments.

| Parameter | Value / Source |
|---|---|
| Grid | 128 × 128 |
| `dt` | 2.5e-5 |
| `tol` | 1e-6 |
| `gam` | `M * (1/128) / (2√2 * atanh(0.9))` |
| `max_it` | `total_time / dt` |
| Output dir | `/project/g_bme-janeslab/SarahG/julia_out/critical_radius_updated_IC/` |
| Job log | `/home/xpz5km/Cahn_Hilliard_Model/Job_specs.csv` |

Initialization: `initialization(nx, nx, method="droplet", h=1/128, R0=R0, gam=gam)` — tanh droplet profile centered in domain.

Output filename:
```
phi_128_<max_it>_1.0e-6__R0_<R0>_eps_<gam>.txt
```

### `julia_run_critical_radius_256.jl`
Same as above but on a 256² grid, allowing finer spatial resolution near the critical radius.

### `julia_run_critical_radius_alt_IC.jl`
Uses an alternate initial condition (e.g., a sharper or differently-centered droplet profile) to test IC sensitivity near the critical radius.

### `julia_run_critical_radius_offcenter.jl`
Places the droplet off-center in the domain to test whether boundary proximity affects critical radius behavior under Neumann BCs.

### `julia_run_critical_radius_periodic.jl`
Runs with periodic boundary conditions using `CH_multigrid_solver_periodic_bc.jl`. Useful to check whether BCs significantly affect the measured critical radius.

### `julia_run_critical_radius_two_drops.jl`
Initializes **two droplets** at x-positions 0.2 and 0.8 on a 512² domain:
```julia
phi = tanh((R0 - R1) / (√2 * gam)) + tanh((R0 - R2) / (√2 * gam)) + 1
```
Guards against R0 ≥ 2.0 (which would cause the droplets to overlap). Output: `/project/g_bme-janeslab/sarah/julia_out/critical_radius_two_drops/`

### `julia_run_critical_radius_two_drops_256_doubledomain.jl`
Two-droplet variant on a 256² grid with a doubled physical domain, giving more separation between droplets.

---

## SLURM Jobs

### `run_julia_multigrid.sh`

| Setting | Value |
|---|---|
| Nodes | 1 |
| Tasks per node | 16 |
| Memory | 50 GB |
| Wall time | 10 hours |
| Array | 1-6 (reads from `options_large_eps_v2.txt`) |
| SLURM output | `./Reports/%A/output.%J.out` |

Each array task reads one line from the options file and passes it as arguments to `julia_run_critical_radius.jl`.

### `run_julia_multigrid_256.sh` and `run_julia_multigrid_alt_IC.sh`
Same structure; differ only in which options file and which `.jl` script they invoke.

---

## Option Files

Each line: `R0  M  total_time`

`options.txt` example:
```
0.100 8 10
0.105 8 10
0.110 8 10
...
```

The files span a range of R0 values bracketing the expected critical radius, with finer sampling near the transition. `M = 8` is the standard interface mesh parameter, giving ε ≈ 0.015 at nx=128.

Versioned files (`v2` through `v9`, `options_large_eps_v2.txt`) represent successive parameter refinements or reruns of failed jobs.

---

## Post-Processing

### `level_set_radius_multiple_droplets.m`
Reads phi output files for multi-droplet runs and tracks the zero level set of each droplet separately, allowing independent radius measurements when the two droplets may dissolve at different times.

### `plot_timepoints_critical_radius.py`
Python script that generates side-by-side heatmap plots of φ at selected timepoints, useful for visually confirming whether a droplet dissolved or grew.

### `plot_timepoints.sh`
Shell launcher for `plot_timepoints_critical_radius.py`.
