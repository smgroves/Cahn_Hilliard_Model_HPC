#!/bin/bash
#SBATCH -N 1
#SBATCH -o ./Reports/compare/%A/output.%j.out 
#SBATCH --ntasks-per-node=16
#SBATCH --account=janeslab
#SBATCH --time=12:00:00
#SBATCH --mem=250G
#SBATCH --partition=standard
#SBATCH --array=1-3

echo $(date)
OPTS=$(sed -n "${SLURM_ARRAY_TASK_ID}"p opts_Julia_redo_v4.txt)
echo $OPTS

# outdir="/project/g_bme-janeslab/SarahG/spinodal_decomp_06_2025/out_MATLAB"
# mkdir -p $outdir
# Load  Julia environment
module load julia/1.9.2
echo "MODULES LOADED"
gridsize=$(echo "$OPTS" | awk '{print $1}')

print_results=$(echo "$OPTS" | awk '{print $3}')
solver=$(echo "$OPTS" | awk '{print $4}')
boundary=$(echo "$OPTS" | awk '{print $2}')
SLURM_ID=${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}

# if [[ "$print_results" == "true" && "$solver" == "NMG" && "$boundary" == "neumann" && "$gridsize" == "512" ]]; then
# if [[ "$print_results" == "true" ]]; then
julia ./rmse.jl $OPTS $SLURM_ID 
# fi


echo "DONE"