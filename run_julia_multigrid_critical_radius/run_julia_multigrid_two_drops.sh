#!/bin/bash
#SBATCH -N 1
#SBATCH -o ./Reports/%A/output.%J.out
#SBATCH --ntasks-per-node=16
#SBATCH --account=janeslab
#SBATCH --time=10:00:00
#SBATCH --mem=50G
#SBATCH --partition=standard
#SBATCH --array=1-32

echo $(date)
OPTS=$(sed -n "${SLURM_ARRAY_TASK_ID}"p options.txt)
outdir="/project/g_bme-janeslab/SarahG/julia_out/critical_radius_two_drops"
mkdir -p $outdir

echo $OPTS
# Load  Julia environment
module load julia/1.9.2
echo "MODULES LOADED"
julia julia_run_critical_radius_two_drops.jl $OPTS

echo "DONE"