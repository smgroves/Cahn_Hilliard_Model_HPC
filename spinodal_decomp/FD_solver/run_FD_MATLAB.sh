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

outdir="/project/g_bme-janeslab/SarahG/spinodal_decomp_04_2025/out_MATLAB/FD"
mkdir -p $outdir


matlab -nodisplay -nosplash -r "run_spinodal_decomp_FD();quit;"
echo "Done."
#3613255