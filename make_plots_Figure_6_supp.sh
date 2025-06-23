#!/bin/bash
#SBATCH -N 1
#SBATCH -o ./Reports/output.%j.out 
#SBATCH --ntasks-per-node=16
#SBATCH --account=janeslab
#SBATCH --time=10:00:00
#SBATCH --mem=50G
#SBATCH --partition=standard


CPC=0.15
cohesin=0.09
nx=512
t=0.04
timesteps=26214
dtout=10
dt=0.000001525878906
frame_rate=100
# # Load  Julia environment
module load julia/1.9.2
module load matlab

echo "MODULES LOADED"

# name="_CPC_0.15_cohesin_0.09_eps_0.0067_smooth"
outdir="/project/g_bme-janeslab/SarahG/julia_out/CPC_geometry/CPC_domain_0_2_e_0.0067_t_0.04"
# mkdir -p $outdir

# julia CPC_geometric_run_array_with_alpha_change_domain_eps_0.0075.jl $CPC $cohesin $nx $t $name


name="phi_${nx}_${timesteps}_1.0e-5__CPC_${CPC}_cohesin_${cohesin}_eps_0.0067_smooth"


# Run Matlab single core program
matlab -nodisplay -r "CHplotting_function_minutes('$outdir', '$name', $dt, $dtout, '_minutes_faster', $frame_rate);quit;"



name="phi_${nx}_${timesteps}_1.0e-5__CPC_${CPC}_cohesin_${cohesin}_eps_0.0067_noisy"

matlab -nodisplay -r "CHplotting_function_minutes('$outdir', '$name', $dt, $dtout, '_minutes_faster', $frame_rate);quit;"

# suffix=""
# matlab -nodisplay -r "plot_single_timepoint_heatmap(1, '$indir','$name', '$suffix');quit;"
# matlab -nodisplay -r "plot_IC_heatmap('$indir','$name', '$suffix', 512, 50);quit;"

echo "Done."
