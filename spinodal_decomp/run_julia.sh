#!/bin/bash
#SBATCH -N 1
#SBATCH -o ./Reports/julia/%A/output.%j.out 
#SBATCH --ntasks-per-node=16
#SBATCH --account=janeslab
#SBATCH --time=8:00:00
#SBATCH --mem=100G
#SBATCH --partition=standard
#SBATCH --array=1

echo $(date)
OPTS=$(sed -n "${SLURM_ARRAY_TASK_ID}"p opts_redo_Julia_v2.txt)
echo $OPTS

outdir="/project/g_bme-janeslab/SarahG/spinodal_decomp_04_2025/out_julia"
mkdir -p $outdir
# Load  Julia environment
module load julia/1.9.2
echo "MODULES LOADED"
# solver=$(echo "$OPTS" | awk '{print $4}')
SLURM_ID=${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}

# if [[ "$solver" == "SAV"  ]]; then
# julia ./CahnHilliard_Julia_solvers/run_spinodal_decomp.jl $OPTS "_25p" $SLURM_ID 
julia ./CahnHilliard_Julia_solvers/run_spinodal_decomp.jl $OPTS $SLURM_ID 

# fi


echo "DONE"