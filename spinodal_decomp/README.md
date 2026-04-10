# spinodal_decomp/

Cross-language reference implementations of the Cahn-Hilliard solver and scripts to compare their output. The primary goal is to verify that Julia (NMG and SAV), MATLAB, and Python solvers produce equivalent results across grid resolutions and boundary conditions.

---

## Directory Structure

```
spinodal_decomp/
├── CahnHilliard_Julia_solvers/     # Julia NMG and SAV solvers (modular, HPC-ready)
├── CahnHilliard_MATLAB_solvers/    # MATLAB solver implementations
├── CahnHilliard_Python_solvers/    # Python solver implementations
├── FD_solver/                      # Simple finite-difference reference solver
├── IC/                             # Shared initial condition files
├── Reports/                        # Simulation output from this study
├── compare_phi_norms.py            # Compute and compare L2 norms of φ across languages
├── compare_phi_mass.py             # Compute and compare mass (integral of φ) across languages
├── compare_L2_error.py             # Compare L2 errors relative to a reference solution
├── compare_norms.csv               # Aggregated norm comparison results
├── compare_norms_v3.csv            # Updated norm comparison (v3)
├── compare_masses.csv              # Aggregated mass comparison results
├── L2_error_Julia.csv              # L2 errors for Julia runs
├── L2_error_MATLAB.csv             # L2 errors for MATLAB runs
├── L2_error_python.csv             # L2 errors for Python runs
├── difference_map_movie.jl         # Julia: generate difference-map movies between two phi fields
├── rmse.jl                         # Julia: compute RMSE between phi fields
├── opts.txt / opts_*.txt           # Parameter option files for batch HPC runs
└── failed_python_06_2025.txt       # Log of failed Python jobs (for resubmission)
```

---

## Julia Solvers (`CahnHilliard_Julia_solvers/`)

Two complete solver implementations, each split into a main driver and helper modules:

### NMG — Nonlinear Multigrid

| File | Purpose |
|---|---|
| `CahnHilliard_NMG.jl` | Main driver; top-level `cahn_hilliard_nmg()` function. Supports Neumann and periodic BCs, outputs phi arrays or text files. |
| `nmg_solver.jl` | V-cycle, restriction, prolongation, and smoothing routines |
| `relax.jl` | Damped Jacobi relaxation for the coupled (φ, μ) block system |
| `laplace.jl` | 5-point Laplacian with switchable BC |
| `ch_error2.jl` | Residual norm computation |
| `ch_initialization.jl` | Initial conditions: random, droplet (tanh), spinodal, CPC geometries |
| `aux_functions_NMG.jl` | Utility functions (mass, energy, output formatting) |

### SAV — Scalar Auxiliary Variable

| File | Purpose |
|---|---|
| `CahnHilliard_SAV.jl` | Main driver; unconditionally energy-stable scheme using a scalar auxiliary variable. Uses FFT-based Laplacian for periodic domains. |
| `sav_solver.jl` | SAV time-stepping: computes the auxiliary variable `r(t)` and updates φ each step |
| `aux_functions_SAV.jl` | SAV-specific utilities (energy computation, auxiliary variable initialization) |

### Shared Utilities

| File | Purpose |
|---|---|
| `run_spinodal_decomp.jl` | Local (non-HPC) run script for spinodal decomposition |
| `run_spinodal_decomp_HPC.jl` | SLURM-ready run script; reads parameters from opts file via `ARGS` |
| `ch_movie.jl` | Generate in-memory movie from phi time series |
| `ch_movie_from_file.jl` | Generate movie by reading saved phi text files |
| `make_movie.m` | MATLAB helper to render saved frames into AVI |

---

## Python Solvers (`CahnHilliard_Python_solvers/`)

Python equivalents of the NMG and/or SAV solvers, used to validate portability on HPC systems where Julia is unavailable.

---

## MATLAB Solvers (`CahnHilliard_MATLAB_solvers/`)

MATLAB implementations of the Cahn-Hilliard solver. Primarily used as the reference solution against which Julia and Python outputs are compared.

---

## FD Solver (`FD_solver/`)

Simple explicit or semi-implicit finite difference solver. Used as a baseline reference, not for production HPC runs.

---

## Cross-Language Comparison Scripts

### `compare_phi_norms.py`
- Walks the `Reports/` directory and finds all `*_phi.csv` files
- Parses solver type, language, grid size, BC, and time step from filenames using regex
- Reads the final timepoint and computes L2 norm of φ
- Outputs aggregated results to `compare_norms.csv`

### `compare_phi_mass.py`
- Same structure as `compare_phi_norms.py` but computes total mass (sum of φ) at each saved timepoint
- Checks mass conservation across languages and methods
- Results saved to `compare_masses.csv`

### `compare_L2_error.py`
- Computes the L2 error of each solver's φ field relative to a designated reference solution (typically high-resolution MATLAB NMG)
- Results saved to `L2_error_*.csv` per language

### `rmse.jl`
- Julia script computing RMSE between two phi fields (e.g., comparing two solver outputs or two timepoints)

### `difference_map_movie.jl`
- Reads two phi output files and generates a frame-by-frame difference map movie to visually diagnose divergence between solvers

---

## Option Files

`opts.txt` and its versioned variants (`opts_06_2025.txt`, `opts_redo_Julia.txt`, etc.) are parameter list files used by HPC array jobs. Each line specifies a set of arguments (grid size, dt, BC, solver type, etc.) that is read by `run_spinodal_decomp_HPC.jl` via `SLURM_ARRAY_TASK_ID`.

The `failed_python_06_2025.txt` and `missing05_04.txt` files log job IDs or parameter sets that did not complete successfully, used when resubmitting partial arrays.

---

## Output Location

Simulation results are written to `Reports/` (local) or to the HPC scratch path specified in the run script (e.g., `/project/g_bme-janeslab/SarahG/julia_out/spinodal_decomp/`). Output filenames encode all run parameters, e.g.:

```
NMG_Julia_2000_dt_2.50e-05_Nx_128_neumann_n_relax_2_phi.csv
```
