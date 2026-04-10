# spinodal_decomp/

Cross-language reference implementations of the Cahn-Hilliard solver and scripts to compare their output. The primary goals are to verify that Julia (NMG and SAV), MATLAB, and Python solvers produce equivalent results across grid resolutions and boundary conditions, and to validate the solvers against the expected t^(1/3) coarsening law via structure factor analysis.

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
├── structure_factor_analysis/      # Structure factor computation and coarsening law validation
├── compare_phi_norms.py            # Compute and compare L2 norms of φ across languages
├── compare_phi_mass.py             # Compute and compare mass (integral of φ) across languages
├── compare_L2_error.py             # Compare L2 errors relative to a reference solution
├── difference_map_movie.jl         # Julia: generate difference-map movies between two phi fields
├── rmse.jl                         # Julia: compute RMSE between Julia phi fields
├── rmse_python.jl                  # Julia: compute RMSE between Python and Julia phi fields
├── compare_norms.csv               # Aggregated norm comparison results
├── compare_norms_v3.csv            # Updated norm comparison (v3)
├── compare_masses.csv              # Aggregated mass comparison results
├── L2_error_Julia.csv              # L2 errors for Julia runs
├── L2_error_MATLAB.csv             # L2 errors for MATLAB runs
├── L2_error_python.csv             # L2 errors for Python runs
├── opts.txt / opts_*.txt           # Parameter option files for batch HPC runs
├── run_julia.sh                    # SLURM: run Julia NMG/SAV solvers
├── run_python.sh                   # SLURM: run Python SAV solvers
├── run_compare_phi_norms.sh        # SLURM: run compare_phi_norms.py
├── run_compare_phi_mass.sh         # SLURM: run compare_phi_mass.py
├── run_compare_L2_error.sh         # SLURM: run compare_L2_error.py
├── run_difference_map_movie.sh     # SLURM array: run difference_map_movie.jl
├── run_rmse.sh                     # SLURM array: run rmse.jl
├── run_rmse_python.sh              # SLURM array: run rmse_python.jl
├── run_sf_main.sh                  # SLURM array: run structure factor analysis
├── monitor_clock_rate_run_julia.sh # SLURM: run Julia solver + log CPU clock rate
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

Python equivalents of the SAV solver, used to validate portability on HPC systems. Run via `run_python.sh` which calls `run_spinodal_decomp_HPC.py` with grid size, boundary condition, solver type, and a SLURM job ID. Output goes to `/project/g_bme-janeslab/SarahG/spinodal_decomp_06_2025/out_python/`.

---

## MATLAB Solvers (`CahnHilliard_MATLAB_solvers/`)

MATLAB implementations of the Cahn-Hilliard solver. Primarily used as the reference solution against which Julia and Python outputs are compared.

---

## FD Solver (`FD_solver/`)

Simple explicit or semi-implicit finite difference solver. Used as a baseline reference, not for production HPC runs.

---

## Structure Factor Analysis (`structure_factor_analysis/`)

Validates that the solvers reproduce the expected Lifshitz-Slyozov-Wagner (LSV) coarsening law l(t) ∝ t^(1/3) by computing the radially-averaged structure factor S(k,t) from the φ field time series.

### `structure_factor.py` — Core Library

**`StructureFactorAnalyzer` class:**

| Method | Description |
|---|---|
| `compute_structure_factor(phi)` | Computes S(k) for a single φ snapshot: mean-subtracts φ, applies 2D FFT, radially averages the power spectrum in pixel-space annuli, converts pixel radius to physical wavenumber k = 2π·r/L. Returns k array, S(k) array, peak wavenumber k_max, and characteristic length l = 2π/k_max. |
| `analyze_series(phi_list, times)` | Runs `compute_structure_factor` over a time series; returns dict with t, k, S(k) list, l(t), and k_peak(t). |

**Standalone functions:**

| Function | Description |
|---|---|
| `fit_exponent(t, l, t_min, t_max)` | Fits l(t) ∝ t^n via log-log linear regression (`np.polyfit`) over an optional time window. Returns exponent, amplitude, R², and fit arrays. |
| `plot_structure_factor_snapshots(result)` | Plots S(k) on a semilog scale at N evenly-spaced timepoints. |
| `plot_coarsening_law(t, l, ...)` | Log-log plot of l(t) with power-law fit (coral dashed) and theory t^(1/3) curve (green dotted). Annotates measured exponent, theory value, and % deviation. |
| `plot_phase_field_samples(phi_list, times)` | Shows φ heatmaps at t=0, midpoint, and final time using a red-blue colormap. |

### `sf_main.py` — Analysis Driver

Accepts command-line arguments: `N  boundary  print_results  method`

For each job:
1. Constructs input filename using the pattern `{method}_Julia_2000_dt_5.50e-06_Nx_{N}_{boundary}_dtout_1050p_phi.csv`
2. Loads phi from `/project/g_bme-janeslab/Kevin_PaperBackups/GrovesSM_PNAS/Data_and_output/simulation_results/spinodal_decomp_06_2025/out_Julia/`
3. Reshapes the flat file into (N, N, T) time series
4. Runs `StructureFactorAnalyzer.analyze_series()` over all T snapshots (times spaced 10 time units apart)
5. Fits the coarsening exponent over the middle 80% of the time range (avoids transient and saturation)
6. Saves three plots and a CSV to `./structure_factor_results/{N}_{boundary}_{method}_{dt}/`:
   - `s_k_evolution.png` — S(k) at 4 snapshots
   - `coarsening_law.png` — l(t) log-log with fit
   - `phase_field_evolution.png` — φ heatmaps at 3 times
   - `fit_results.csv` — exponent, amplitude, R², fit window

Prints a pass/fail summary: ✓ if |exponent - 1/3| < 0.01, ≈ if < 0.05, ✗ otherwise.

### `run_sf_main.sh` — SLURM Job

| Setting | Value |
|---|---|
| Memory | 500 GB |
| Wall time | 30 hours |
| Array size | 1–32 (reads from `opts.txt`) |
| SLURM output | `./Reports/python/%A/output.%j.out` |
| Module | `apptainer`, `miniforge` |

Option file format: `N  boundary  print  solver`

---

## Cross-Language Comparison Scripts

### `compare_phi_norms.py`
Walks the `Reports/` directory, finds all `*_phi.csv` files, parses solver/language/grid/BC from filenames, reads the final timepoint, computes L2 norm of φ, and outputs `compare_norms.csv`.

### `compare_phi_mass.py`
Same structure; computes total mass (sum of φ) at each saved timepoint to check mass conservation. Results in `compare_masses.csv`.

### `compare_L2_error.py`
Computes the L2 error of each solver's φ field relative to a designated reference solution. Results in `L2_error_*.csv` per language.

### `rmse.jl`
Computes RMSE between two Julia phi fields (e.g., NMG vs. SAV, or different grid sizes). Reads parameters from `opts_Julia_redo_v4.txt` via `ARGS`.

### `rmse_python.jl`
Computes RMSE between Python SAV and Julia SAV phi fields for validation. Reads from `opts_06_2025.txt`.

### `difference_map_movie.jl`
Generates a frame-by-frame difference map movie between two phi output files (e.g., Julia vs. Python SAV) to visually diagnose where and when solvers diverge.

---

## SLURM Scripts Summary

| Script | What it runs | Array size | Memory |
|---|---|---|---|
| `run_julia.sh` | `run_spinodal_decomp_HPC.jl` | 1–96 | 1400 GB |
| `run_python.sh` | `run_spinodal_decomp_HPC.py` (SAV only) | 1–96 | 100 GB |
| `run_compare_phi_norms.sh` | `compare_phi_norms.py` | single job | 250 GB |
| `run_compare_phi_mass.sh` | `compare_phi_mass.py` | single job | 250 GB |
| `run_compare_L2_error.sh` | `compare_L2_error.py` | single job | 250 GB |
| `run_difference_map_movie.sh` | `difference_map_movie.jl` (SAV only) | 1–96 | 250 GB |
| `run_rmse.sh` | `rmse.jl` | 1–3 | 250 GB |
| `run_rmse_python.sh` | `rmse_python.jl` (SAV only) | 1–96 | 250 GB |
| `run_sf_main.sh` | `sf_main.py` | 1–32 | 500 GB |
| `monitor_clock_rate_run_julia.sh` | `run_spinodal_decomp_HPC.jl` + CPU freq monitor | 1–3 | 500 GB |

`monitor_clock_rate_run_julia.sh` runs a background loop logging `grep MHz /proc/cpuinfo` every 2 seconds alongside the Julia simulation. Saves per-job CPU frequency logs to `out_julia/cpu_freq_<SLURM_ID>.log` and prints the average MHz at the end.

---

## Option Files

`opts.txt` and versioned variants (`opts_06_2025.txt`, `opts_Julia_redo_v4.txt`, etc.) are parameter list files. Each line: `N  boundary  print  solver`

Example:
```
128 neumann true NMG
256 periodic true SAV
512 neumann true NMG
```

`failed_python_06_2025.txt` and `missing05_04.txt` log job IDs or parameter sets that did not complete, used when resubmitting partial arrays.

---

## Output Locations

| Data | Path |
|---|---|
| Julia phi output | `/project/g_bme-janeslab/SarahG/spinodal_decomp_06_2025/out_julia/` |
| Python phi output | `/project/g_bme-janeslab/SarahG/spinodal_decomp_06_2025/out_python/` |
| Paper backup (Julia) | `/project/g_bme-janeslab/Kevin_PaperBackups/GrovesSM_PNAS/Data_and_output/simulation_results/spinodal_decomp_06_2025/out_Julia/` |
| Structure factor results | `./structure_factor_results/{N}_{boundary}_{method}_{dt}/` |
| SLURM logs | `./Reports/{julia,python,compare}/` |

Output phi filenames encode all run parameters, e.g.:
```
NMG_Julia_2000_dt_5.50e-06_Nx_512_neumann_dtout_1050p_phi.csv
```
