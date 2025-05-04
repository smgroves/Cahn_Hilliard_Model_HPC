#!/bin/bash
#SBATCH -N 1
#SBATCH -o ./Reports/%A/output.%J.out
#SBATCH --ntasks-per-node=16
#SBATCH --account=janeslab
#SBATCH --time=10:00:00
#SBATCH --mem=50G
#SBATCH --partition=standard
#SBATCH --array=1-66

echo $(date)
OPTS=$(sed -n "${SLURM_ARRAY_TASK_ID}"p options_alpha_v6.txt)
echo $OPTS
module load matlab
echo "LOADED MATLAB"
R0=$(echo "$OPTS" | awk '{print $1}')
m=$(echo "$OPTS" | awk '{print $2}')
alpha=$(echo "$OPTS" | awk '{print $4}')

echo $R0
echo $m
echo $alpha

matlab -nodisplay -r "level_set_radius_array_alpha($R0, $m, $alpha); exit;"
echo "DONE"
