#!/bin/bash
#SBATCH -N 1
#SBATCH -o ./Reports/%A/output.%J.out
#SBATCH --ntasks-per-node=16
#SBATCH --account=janeslab
#SBATCH --time=30:00:00
#SBATCH --mem=200G
#SBATCH --partition=standard
#SBATCH --array=1-5

echo $(date)
OPTS=$(sed -n "${SLURM_ARRAY_TASK_ID}"p options_periodic_check.txt)
echo $OPTS
# Load  Julia environment
module load julia/1.9.2
echo "MODULES LOADED"
julia julia_run_critical_radius_256.jl $OPTS

echo "DONE"