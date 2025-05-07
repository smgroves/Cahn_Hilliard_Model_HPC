#!/bin/bash
#SBATCH -N 1
#SBATCH -o ./Reports/output.%j.out 
#SBATCH --ntasks-per-node=16
#SBATCH --account=janeslab
#SBATCH --time=10:00:00
#SBATCH --mem=250G
#SBATCH --partition=standard

# indir="/scratch/xpz5km/Cahn_Hilliard_Model/julia_out/critical_radius/CPC_geometry"
# indir="/project/g_bme-janeslab/SarahG/julia_out/CPC_geometry/CPC_domain_0_2_e_0.0075_noisy_cohesin/sd_0.11/individual_seeds"
# indir="/project/g_bme-janeslab/SarahG/julia_out/critical_radius_alt_IC"

indir="/project/g_bme-janeslab/SarahG/julia_out/CPC_geometry/CPC_domain_0_2_e_0.0067_t_0.4"

dtout=10
dt=0.000001525878906
nx=512
frame_rate=20
# total_time=.03
# dtdec=$(printf "%.14f" $dt)
# timesteps=$(echo " $total_time / $dtdec" | bc)
# timesteps=163840
# for cohesin in 2 4
# do

# name="phi_512_19661_1.0e-5__CPC_0.22_cohesin_0.09_eps_0.0075_alpha_0_domain_0_2"
# name="phi_128_400000_1.0e-6__R0_0.149_eps_0.030019_two_halves"
name="phi_512_262144_1.0e-5__CPC_0.25_cohesin_0.11_eps_0.0067_domain_0_2"

echo $(date)
module load matlab
echo "LOADED MATLAB"

# Run Matlab single core program
matlab -nodisplay -r "CHplotting_function('$indir', '$name', $dt, $dtout, '', $frame_rate);quit;"

# suffix=""
# matlab -nodisplay -r "plot_single_timepoint_heatmap(1, '$indir','$name', '$suffix');quit;"
# matlab -nodisplay -r "plot_IC_heatmap('$indir','$name', '$suffix', 512, 50);quit;"

echo "Done."
