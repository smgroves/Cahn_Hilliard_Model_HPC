#!/bin/bash
#SBATCH -N 1
#SBATCH -o ../Reports/%A/output.%j.out 
#SBATCH --ntasks-per-node=16
#SBATCH --account=janeslab
#SBATCH --time=10:00:00
#SBATCH --mem=300G
#SBATCH --partition=standard
#SBATCH --array=1-10

export MATLAB_DISABLE_SERVICEHOST=true
module load matlab
# seed="1111"
# indir="/project/g_bme-janeslab/SarahG/julia_out/CPC_geometry/CPC_domain_0_2_e_0.0075_noisy_cohesin/sd_0.25/seed_${seed}"
# indir="/project/g_bme-janeslab/SarahG/julia_out/CPC_geometry/CPC_domain_0_2_e_0.0075_alpha_-0.5_new_IC/"
indir="/project/g_bme-janeslab/SarahG/julia_out/CPC_geometry/CPC_domain_0_2_e_0.0067_noisy_cohesin/sd_0.11/"
# indir="/project/g_bme-janeslab/SarahG/julia_out/CPC_geometry/CPC_domain_0_2_e_0.0089_noisy_cohesin/"

# OPTS=$(sed -n "${SLURM_ARRAY_TASK_ID}"p CPC_geometric_array_with_alpha_change_domain_options_eps_0.0075_v7.txt)
# indir="/project/g_bme-janeslab/SarahG/julia_out/CPC_geometry/CPC_domain_0_2_e_0.0067"
# indir="/project/g_bme-janeslab/SarahG/julia_out/CPC_geometry/CPC_domain_0_2_e_0.0067_t_0.04"
# indir="/project/g_bme-janeslab/SarahG/julia_out/CPC_geometry/CPC_domain_0_2_e_0.0067_t_0.4"

OPTS=$(sed -n "${SLURM_ARRAY_TASK_ID}"p CPC_geometric_array_with_alpha_change_domain_options_eps_0.0075_v7.txt)
# OPTS=$(sed -n "${SLURM_ARRAY_TASK_ID}"p CPC_geometric_array_with_alpha_change_domain_options_eps_0.0067_large_widths_t_0.4.txt)

# OPTS=$(sed -n "${SLURM_ARRAY_TASK_ID}"p CPC_geometric_array_with_alpha_change_domain_options_eps_0.0075_alpha_-0.5_v2.txt)
dtout=10
dt=0.000001525878906
nx=512
CPC=$(echo "$OPTS" | awk '{print $1}')
cohesin=$(echo "$OPTS" | awk '{print $2}')
time=$(echo "$OPTS" | awk '{print $4}')

echo $CPC
echo $cohesin
echo $time
# Check if the entered number is even
if [ $time == 0.015 ]; then
timesteps=9830
elif [ $time == 0.03 ]; then
timesteps=19661
elif [ $time == 0.045 ]; then
timesteps=29491
elif [ $time == 0.4 ]; then
timesteps=262144
fi

echo $timesteps
# if [ $timesteps == 29491 ]; then
name="phi_${nx}_${timesteps}_1.0e-5__CPC_${CPC}_cohesin_${cohesin}_eps_0.0067_alpha_0_domain_0_2"
# name = sprintf('phi_%s_%s_1.0e-5__CPC_%s_cohesin_%s_eps_0.0067_alpha_0_domain_0_2',string(nx),string(timesteps),CPC,cohesin)
echo $name 
outdir="$indir/$name"
echo $outdir
# outdir=sprintf('%s%s',indir,name)
# if [[ "$CPC" == "0.35" || "$CPC" == "0.2" ]]; then
# #2112907
#     echo "CPC is either 0.35 or 0.2"
matlab -nodisplay -nosplash -r "kymograph_central_droplets_domain('$indir', '$outdir','$name', $dt, $dtout, false, 0, 6.4);quit;"
    # echo "Done."
# fi

