#!/bin/bash
#SBATCH -N 1
#SBATCH -o ../Reports/MATLAB/%A/output.%j.out 
#SBATCH --ntasks-per-node=16
#SBATCH --account=janeslab
#SBATCH --time=25:00:00
#SBATCH --mem=100G
#SBATCH --partition=standard
#SBATCH --array=1

echo $(date)
module load matlab

OPTS=$(sed -n "${SLURM_ARRAY_TASK_ID}"p ../opts_redo_MATLAB_v2.txt)
echo $OPTS

outdir="/project/g_bme-janeslab/SarahG/spinodal_decomp_04_2025/out_MATLAB"
mkdir -p $outdir
# Load  Julia environment
module load julia/1.9.2
echo "MODULES LOADED"
GridSize=$(echo "$OPTS" | awk '{print $1}')
boundary=$(echo "$OPTS" | awk '{print $2}')
print_results=$(echo "$OPTS" | awk '{print $3}')
solver=$(echo "$OPTS" | awk '{print $4}')

SLURM_ID=${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}

matlab -nodisplay -nosplash -r "run_spinodal_decomp($GridSize, '$boundary', '$print_results', '$solver', '$SLURM_ID', '_25p');quit;"
echo "Done."