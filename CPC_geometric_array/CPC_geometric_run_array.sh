#!/bin/bash
#SBATCH -N 1
#SBATCH -o ../Reports/output.%j.out 
#SBATCH --ntasks-per-node=16
#SBATCH --account=janeslab
#SBATCH --time=10:00:00
#SBATCH --mem=50G
#SBATCH --partition=standard
#SBATCH --array=1-30

echo $(date)
OPTS=$(sed -n "${SLURM_ARRAY_TASK_ID}"p CPC_options_eps_0.039117.txt)
echo $OPTS
# Load  Julia environment
module load julia/1.9.2
echo "MODULES LOADED"
julia CPC_geometric_run_array.jl $OPTS

echo "DONE"