#!/bin/bash
#SBATCH -N 1
#SBATCH -o ./Reports/compare/output.%J.out
#SBATCH --ntasks-per-node=16
#SBATCH --account=janeslab
#SBATCH --time=10:00:00
#SBATCH --mem=250G
#SBATCH --partition=standard

#3035034
#4156068
echo $(date)
# Load  Julia environment
module load miniforge/24.3.0-py3.11
echo "MODULES LOADED"
# out_file=f"/home/xpz5km/Cahn_Hilliard_Model_HPC/spinodal_decomp/compare_norms.csv"

python compare_L2_error.py

echo "DONE"