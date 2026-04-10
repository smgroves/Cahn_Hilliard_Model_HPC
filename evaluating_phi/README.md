# evaluating_phi/

Scripts for analyzing how the order parameter φ behaves under different geometries and free-energy asymmetry values (α). Focuses on characterizing equilibrium states, maximum/minimum φ values, and interface profiles for both flat-interface and single-droplet configurations.

---

## Directory Structure

```
evaluating_phi/
├── Reports/                              # Simulation output files used as input here
├── f_with_alpha_range.m                  # MATLAB: max/min φ vs. time for CPC geometry, α = -0.5
├── f_with_alpha_range_single_droplet.m   # MATLAB: max/min φ vs. time for single droplet, varying α
├── f_with_alpha_flat_interface.m         # MATLAB: equilibrium φ profile for a flat interface
├── flat_interface_eq_interface_sigmoid.m # MATLAB: sigmoid fit to flat interface equilibrium profile
├── run_f_with_alpha.sh                   # Shell: MATLAB batch launcher for alpha range scripts
└── run_eq_interface_sigmoid.sh           # Shell: MATLAB batch launcher for sigmoid fit script
```

---

## Scripts

### `f_with_alpha_range.m`
Reads CPC geometry phi output files from the Reports directory (path: `/project/g_bme-janeslab/sarah/alpha_nonzero/CPC_geometry/CPC_alpha_-0.5`) for α = -0.5. For each combination of `CPC_width` and `cohesin_width`, it:
1. Loads and reshapes the concatenated phi text file into a 3D array (x, y, time)
2. Computes `max(φ)` and `min(φ)` at each saved timepoint
3. Plots these vs. time on two separate figures
4. Saves output PDFs to the input directory

The equilibrium max/min φ values are used to verify that the modified double-well potential (with α ≠ 0) produces the expected asymmetric phase equilibria.

### `f_with_alpha_range_single_droplet.m`
Same approach as above but for single-droplet initial conditions over a range of α values. Useful for isolating the effect of α on bulk phase values independent of geometry.

### `f_with_alpha_flat_interface.m`
Analyzes phi profiles for flat interface simulations to extract the equilibrium interface shape. Used to validate that the solver reproduces the expected tanh profile at steady state.

### `flat_interface_eq_interface_sigmoid.m`
Fits a sigmoid (logistic) function to the equilibrium flat-interface φ profile. The fit parameters (midpoint, width) can be used to back out the effective interface thickness and compare it against the theoretical value set by ε.

---

## Shell Launchers

`run_f_with_alpha.sh` and `run_eq_interface_sigmoid.sh` call MATLAB in batch mode (`matlab -nodisplay -nosplash -r`) to run the corresponding `.m` scripts on the HPC cluster without a GUI.

---

## Input Data

All scripts read phi output files from `Reports/` or from the HPC project storage path embedded in the script. Files follow the naming convention:

```
phi_256_2000_1.0e-5__CPC_<width>_cohesin_<width>_eps_0.04_alpha_-0.5.txt
```

Each file contains the φ field at all saved timepoints concatenated row-wise; scripts reshape using `phidims` computed from the file dimensions.

---

## Output

Figures are saved as PDFs directly to the input data directory. No new simulation output is generated here.
