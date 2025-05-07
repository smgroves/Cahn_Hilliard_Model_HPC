#!/bin/bash
#SBATCH -N 1
#SBATCH -o ./Reports/output.%j.out 
#SBATCH --ntasks-per-node=16
#SBATCH --account=janeslab
#SBATCH --time=26:00:00
#SBATCH --mem=400G
#SBATCH --partition=standard

export MATLAB_DISABLE_SERVICEHOST=true
dtout=10
dt=0.000001525878906

CPC=0.15
cohesin=0.09
nx=512
t=0.04
timesteps=26214

# # Load  Julia environment
module load julia/1.9.2
module load matlab

echo "MODULES LOADED"

# name="_CPC_0.15_cohesin_0.09_eps_0.0067_smooth"
outdir="/project/g_bme-janeslab/SarahG/julia_out/CPC_geometry/CPC_domain_0_2_e_0.0067_t_0.04"
# mkdir -p $outdir

# julia CPC_geometric_run_array_with_alpha_change_domain_eps_0.0075.jl $CPC $cohesin $nx $t $name


name="phi_${nx}_${timesteps}_1.0e-5__CPC_${CPC}_cohesin_${cohesin}_eps_0.0067_smooth"
outdir_plot="$outdir/$name"
echo $outdir_plot

matlab -nodisplay -nosplash -r "kymograph_central_droplets_domain('$outdir', '$outdir_plot','$name', $dt, $dtout, false, 0, 6.4);quit;"


# name="_CPC_0.15_cohesin_0.09_eps_0.0067_noisy"
# julia CPC_geometric_run_array_with_alpha_change_domain_eps_0.0075_noisy_cohesin.jl $CPC $cohesin $nx $t 42 $name

name="phi_${nx}_${timesteps}_1.0e-5__CPC_${CPC}_cohesin_${cohesin}_eps_0.0067_noisy"
outdir_plot="$outdir/$name"
echo $outdir_plot

matlab -nodisplay -nosplash -r "kymograph_central_droplets_domain('$outdir', '$outdir_plot','$name', $dt, $dtout, false, 0, 6.4);quit;"


echo "DONE"