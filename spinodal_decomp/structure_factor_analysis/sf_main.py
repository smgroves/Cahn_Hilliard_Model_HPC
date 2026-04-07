import matplotlib.font_manager as fm
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from matplotlib import font_manager
from matplotlib import rcParams
 
import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft2, fftfreq
from scipy.ndimage import gaussian_filter
from scipy.optimize import curve_fit
import sys
import os
 
# Import structure_factor module (adjust path if needed)
sys.path.insert(0, './structure_factor_analysis')
from structure_factor import *
 

N = int(sys.argv[1])
boundary = sys.argv[2]
print_results = sys.argv[3] == "true"
method = sys.argv[4]

dt = "5.50e-06"

outdir=f"./structure_factor_results/{N}_{boundary}_{method}_{dt}"
os.mkdir(outdir)

plt.rcParams['pdf.use14corefonts'] = True
 
 
print(f"Running analysis: N={N}, boundary={boundary}, method={method}")
 
L = 1
 
# Construct input filename based on parameters
# Adjust path and naming convention to match your actual files
indir = f"/project/g_bme-janeslab/Kevin_PaperBackups/GrovesSM_PNAS/Data_and_output/simulation_results/spinodal_decomp_06_2025/out_Julia"
phi_filename = f"{method}_Julia_2000_dt_{dt}_Nx_{N}_{boundary}_dtout_1050p_phi.csv"
phi_path = os.path.join(indir, phi_filename)
#  SAV_Julia_20000_dt_5.50e-07_Nx_512_neumann_dtout_1075p_energy
print(f"Loading: {phi_path}")
if not os.path.exists(phi_path):
    print(f"ERROR: File not found: {phi_path}")
    sys.exit(1)
 
phi_raw = np.genfromtxt(phi_path)
phi_raw = phi_raw.reshape(-1, N, N).transpose(1, 2, 0)
 
print(f"Loaded shape: {phi_raw.shape}")
 
# %%
phi_list, times = [], []
 
for i in range(phi_raw.shape[2]):
    phi_list.append(phi_raw[:, :, i])
times = np.arange(phi_raw.shape[2]) * 10.0
 
print(
    f"Generated {len(phi_list)} snapshots over t ∈ [{times[0]:.2f}, {times[-1]:.2f}]"
)
 
# %%
# Analyze
print("Computing structure factor S(k,t)...")
analyzer = StructureFactorAnalyzer(L=L, N=N)
result = analyzer.analyze_series(phi_list, times=times)
print(
    f"Length scales range: l(t) ∈ [{result['l'].min():.3f}, {result['l'].max():.3f}]"
)
 
# %%
# Fit
print("\n[3] Fitting coarsening exponent...")
# Use middle 80% of data to avoid transient and saturation
t_min = times[int(0.2 * len(times))]
t_max = times[-1]
 
fit = fit_exponent(result['t'], result['l'], t_min=t_min, t_max=t_max)
 
print(f"Exponent (fit):    {fit['exponent']:.6f}")
print(f"Exponent (theory): 0.333333")
print(
    f"Error:             {abs(fit['exponent'] - 1/3) / (1/3) * 100:.2f}%"
)
print(f"Amplitude:         {fit['amplitude']:.6f}")
print(f"R2 quality:        {fit['r2']:.8f}")
print(f"Fit window:        t ∈ [{t_min:.3f}, {t_max:.3f}]")
 
# %%
# Plots
print("\n[4] Generating plots...")
 
fig1, _ = plot_structure_factor_snapshots(result, n_snapshots=4)
fig1.savefig(os.path.join(outdir, 's_k_evolution.png'),
             dpi=150,
             bbox_inches='tight')
print(f"✓ Saved: {outdir}/s_k_evolution.png")
 
fig2, _, _ = plot_coarsening_law(result['t'],
                                 result['l'],
                                 t_min=t_min,
                                 t_max=t_max)
fig2.savefig(os.path.join(outdir, 'coarsening_law.png'),
             dpi=150,
             bbox_inches='tight')
print(f"✓ Saved: {outdir}/coarsening_law.png")
 
fig3, _ = plot_phase_field_samples(phi_list, times=times)
fig3.savefig(os.path.join(outdir, 'phase_field_evolution.png'),
             dpi=150,
             bbox_inches='tight')
print(f"✓ Saved: {outdir}/phase_field_evolution.png")
 
# Save fit results to CSV
fit_results = pd.DataFrame({
    'N': [N],
    'boundary': [boundary],
    'method': [method],
    'exponent': [fit['exponent']],
    'exponent_theory': [1/3],
    'error_percent': [abs(fit['exponent'] - 1/3) / (1/3) * 100],
    'amplitude': [fit['amplitude']],
    'r2': [fit['r2']],
    't_min': [t_min],
    't_max': [t_max]
})
fit_results.to_csv(os.path.join(outdir, 'fit_results.csv'), index=False)
print(f"✓ Saved: {outdir}/fit_results.csv")
 
# Summary
print("\n" + "="*60)
if abs(fit['exponent'] - 1 / 3) < 0.01:
    print("✓ SUCCESS: Exponent is consistent with t^(1/3) LSV scaling!")
elif abs(fit['exponent'] - 1 / 3) < 0.05:
    print("≈ MARGINAL: Exponent is close but shows some deviation.")
    print("  Check transient regime and domain saturation effects.")
else:
    print("✗ FAILURE: Exponent significantly deviates from t^(1/3).")
    print("  Likely issues:")
    print("    - Too much numerical dissipation or coarse grid")
    print("    - Data dominated by transient or saturation regime")
    print("    - Non-conservative dynamics (check model type)")
print("="*60)
 
plt.close('all')