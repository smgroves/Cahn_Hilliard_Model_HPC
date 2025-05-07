#!/bin/bash
#SBATCH -N 1
#SBATCH -o ./Reports/%A/output.%j.out 
#SBATCH --ntasks-per-node=16
#SBATCH --account=janeslab
#SBATCH --time=10:00:00
#SBATCH --mem=50G
#SBATCH --partition=standard
#SBATCH --array=1-45

module load matlab
# seed="1111"
# indir="/project/g_bme-janeslab/SarahG/julia_out/CPC_geometry/CPC_domain_0_2_e_0.0075_noisy_cohesin/sd_0.25/seed_${seed}"
# indir="/project/g_bme-janeslab/SarahG/julia_out/CPC_geometry/CPC_domain_0_2_e_0.0075_noisy_cohesin/sd_0.11/individual_seeds"
indir="/project/g_bme-janeslab/SarahG/julia_out/CPC_geometry/CPC_domain_0_1_e_0.0075"

OPTS=$(sed -n "${SLURM_ARRAY_TASK_ID}"p CPC_geometric_array_eps_0.0075.txt)
contour_level=0
dtout=10
dt=0.000001525878906
nx=256
alpha="0"

CPC=$(echo "$OPTS" | awk '{print $1}')
cohesin=$(echo "$OPTS" | awk '{print $2}')
time=$(echo "$OPTS" | awk '{print $4}')
echo $CPC
echo $cohesin
echo $time
matlab -nodisplay -nosplash -r "cd ..; level_set_radius_multiple_droplets($CPC, $cohesin, 0.0075, '$indir', '$alpha', $nx, 0.000001525878906, $time, '_domain_0_1', $contour_level);quit;"
echo "Done."





