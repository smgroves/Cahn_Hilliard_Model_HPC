# CPC_1024_geometry/

High-resolution (1024²) CPC geometry simulation for the manuscript. This folder runs a single production-quality simulation of a round CPC domain at 4× the resolution of the standard 256² runs, and includes a restart capability for extending the run if it hits the wall-time limit.

---

## Directory Structure

```
CPC_1024_geometry/
├── run_solver_manuscript_1024.jl   # Main run script (fresh start)
├── restart_CPC_1024.jl             # Restart run script (resumes from last saved phi)
├── CPC_geometric_run_1024.sh       # SLURM script for fresh run
├── restart_CPC_1024.sh             # SLURM script for restart
└── last_n_lines.py                 # Utility: extract the last N lines of a phi output file
```

---

## Run Parameters (`run_solver_manuscript_1024.jl`)

| Parameter | Value |
|---|---|
| Grid size | 1024 × 1024 |
| Tolerance | 1e-4 |
| Time step `dt` | `0.1 * (1/256)²` ≈ 1.525e-6 |
| Interface thickness `gam` | computed as M=8 at nx=256, ≈ 0.015 |
| Max iterations | 19660 (≈ twice the standard 256² run count for the same physical time) |
| CPC_width | 40 grid points |
| cohesin_width | 16 grid points |
| Output stride `ns` | 10 (save every 10th step) |

The 1024² grid uses the same physical domain (0 to 1) and the same ε as the 256² runs (M=8 at the 256-point scale), so the interface is represented by more grid points at 1024 but models the same physical thickness.

### Initialization
Calls `initialize_round_CPC(1024, 1024, CPC_width=40, cohesin_width=16)` from `../CH_multigrid_solver.jl`. The ring widths are scaled proportionally from the 256² values (CPC=10, cohesin=4 → CPC=40, cohesin=16 at 4×).

### Output
- Phi field written to: `/scratch/xpz5km/Cahn_Hilliard_Model/julia_out/CPC_geometry_1024/`
- Output filename: `phi_1024_19660_1.0e-4__CPC_40_cohesin_16_eps_<gam>.txt`
- Timing appended to: `/home/xpz5km/Cahn_Hilliard_Model/Job_specs.csv`

---

## SLURM Job: `CPC_geometric_run_1024.sh`

| Setting | Value |
|---|---|
| Nodes | 1 |
| Tasks per node | 16 |
| Memory | 100 GB |
| Wall time | 12 hours |
| Email notifications | begin, end, fail → xpz5km@virginia.edu |
| SLURM output | `./Reports/output.%J.out` |

---

## Restart Capability

If the job hits the 12-hour wall-time limit before completing all iterations, `restart_CPC_1024.jl` can resume from the last saved phi frame.

### `last_n_lines.py`
A Python utility that extracts the last N lines from a large phi output file. This is used to recover the final saved phi frame from an incomplete run, which is then used as the initial condition for the restart.

Usage:
```bash
python last_n_lines.py <phi_file> <N> <output_file>
```
where N = nx (e.g., 1024 for a 1024² grid), since each row of the phi file is one row of the spatial grid.

### `restart_CPC_1024.jl`
Reads the last phi frame from the existing output file, uses it as the initial condition, and continues the simulation for the remaining iterations. Appends new phi frames to the same output file.

### `restart_CPC_1024.sh`
SLURM script for the restart job, identical resources to the original run.
