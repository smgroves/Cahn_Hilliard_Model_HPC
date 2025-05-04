#!/bin/bash
#SBATCH -N 1
#SBATCH -o ./Reports/%A/output.%J.out
#SBATCH --ntasks-per-node=16
#SBATCH --account=janeslab
#SBATCH --time=10:00:00
#SBATCH --mem=150G
#SBATCH --partition=standard
#SBATCH --array=1-78

echo $(date)
OPTS=$(sed -n "${SLURM_ARRAY_TASK_ID}"p options_alpha_+0.5_v2.txt)
module load matlab
echo "LOADED MATLAB"
R0=$(echo "$OPTS" | awk '{print $1}')
epsilon=$(echo "$OPTS" | awk '{print $2}')
alpha=$(echo "$OPTS" | awk '{print $4}')
nx=$(echo "$OPTS" | awk '{print $5}')
echo $R0
echo $epsilon
echo $alpha
echo $nx
indir="/project/g_bme-janeslab/SarahG/julia_out/critical_radius_alpha_updated"

matlab -nodisplay -r "level_set_radius_array_alpha_v3($R0, $epsilon, $alpha, $nx, '$indir', 2.5e-5, +0.1666); exit;"
echo "DONE"
