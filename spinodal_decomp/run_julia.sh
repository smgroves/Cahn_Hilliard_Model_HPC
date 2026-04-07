#!/bin/bash
#SBATCH -N 1
#SBATCH -o ./Reports/julia/%A/output.%j.out 
#SBATCH --ntasks-per-node=16
#SBATCH --account=janeslab
#SBATCH --time=8:00:00
#SBATCH --mem=1400G
#SBATCH --partition=standard
#SBATCH --array=1-96

echo $(date)
OPTS=$(sed -n "${SLURM_ARRAY_TASK_ID}"p opts_Julia_redo_v4.txt)
echo $OPTS

outdir="/project/g_bme-janeslab/SarahG/spinodal_decomp_06_2025/out_julia"
mkdir -p $outdir
# Load  Julia environment
module load julia/1.11.6
echo "MODULES LOADED"
boundary=$(echo "$OPTS" | awk '{print $2}')
print=$(echo "$OPTS" | awk '{print $3}')

SLURM_ID=${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}

if [[ "$print" == "true"  ]]; then
# julia ./CahnHilliard_Julia_solvers/run_spinodal_decomp.jl $OPTS "_25p" $SLURM_ID 
    julia ./CahnHilliard_Julia_solvers/run_spinodal_decomp_HPC.jl $OPTS $SLURM_ID 
fi


echo "DONE"