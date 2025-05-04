#!/bin/bash
#SBATCH -N 1
#SBATCH -o ./Reports/%A/output.%j.out 
#SBATCH --ntasks-per-node=16
#SBATCH --account=janeslab
#SBATCH --time=10:00:00
#SBATCH --mem=50G
#SBATCH --partition=standard
#SBATCH --array=1-32

module load matlab
# seed="1111"
# indir="/project/g_bme-janeslab/SarahG/julia_out/
indir="/project/g_bme-janeslab/SarahG/julia_out/critical_radius_two_drops"

OPTS=$(sed -n "${SLURM_ARRAY_TASK_ID}"p options.txt)
contour_level=0
dtout=10
dt=2.5e-5
nx=128
alpha="0"
R0=$(echo "$OPTS" | awk '{print $1}')
m=$(echo "$OPTS" | awk '{print $2}')

if [ $m == 8 ]; then
epsilon="0.015009"
elif [ $m == 16 ]; then
epsilon="0.030019"
elif [ $m == 32 ]; then
epsilon="0.060037"
fi

time=$(echo "$OPTS" | awk '{print $3}')

echo $R0
echo $time
# phi_512_19661_1.0e-5__8_CPC_0.3_cohesin_0.09_eps_0.01_alpha_0_domain_0_2

# phi_128_400000_1.0e-6__R0_0.18_eps_0.060037_two_drops
matlab -nodisplay -nosplash -r "level_set_radius_multiple_droplets($R0, '$epsilon', '$indir', '$alpha', $nx, $dt, $time, '_two_drops', $contour_level);quit;"
echo "Done."





