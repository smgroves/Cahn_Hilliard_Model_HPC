#!/bin/bash
#SBATCH -N 1
#SBATCH -o ./Reports/%A/output.%J.out
#SBATCH --ntasks-per-node=16
#SBATCH --account=janeslab
#SBATCH --time=10:00:00
#SBATCH --mem=200G
#SBATCH --partition=standard
#SBATCH --array=1-32

echo $(date)
OPTS=$(sed -n "${SLURM_ARRAY_TASK_ID}"p ./run_julia_multigrid_critical_radius/options.txt)
echo $OPTS
module load matlab
echo "LOADED MATLAB"
R0=$(echo "$OPTS" | awk '{print $1}')
m=$(echo "$OPTS" | awk '{print $2}')
echo $R0
echo $m
suffix="_periodic"
# indir="/scratch/xpz5km/Cahn_Hilliard_Model/julia_out/critical_radius_updated_IC"
indir="/project/g_bme-janeslab/SarahG/julia_out/critical_radius_periodic"
matlab -nodisplay -r "level_set_radius_array_alt_IC($R0, $m, 128,'$indir', '$suffix'); exit;"
echo "DONE"
