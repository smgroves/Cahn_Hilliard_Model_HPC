#!/bin/bash
#SBATCH -N 1
#SBATCH -o ./Reports/compare/%A/output.%j.out 
#SBATCH --ntasks-per-node=16
#SBATCH --account=janeslab
#SBATCH --time=12:00:00
#SBATCH --mem=250G
#SBATCH --partition=standard
#SBATCH --array=1-96

echo $(date)
OPTS=$(sed -n "${SLURM_ARRAY_TASK_ID}"p opts_06_2025.txt)
echo $OPTS

# outdir="/project/g_bme-janeslab/SarahG/spinodal_decomp_06_2025/out_MATLAB"
# mkdir -p $outdir
# Load  Julia environment
module load julia/1.9.2
echo "MODULES LOADED"
print_results=$(echo "$OPTS" | awk '{print $3}')
solver=$(echo "$OPTS" | awk '{print $4}')
boundary=$(echo "$OPTS" | awk '{print $2}')
SLURM_ID=${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}

if [[ "$print_results" == "true" && "$solver" == "SAV" ]]; then
    julia ./rmse.jl $OPTS $SLURM_ID 
fi


echo "DONE"