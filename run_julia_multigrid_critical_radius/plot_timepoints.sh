#!/bin/bash
#SBATCH -N 1
#SBATCH -o ./Reports/%A/output.%J.out
#SBATCH --ntasks-per-node=16
#SBATCH --account=janeslab
#SBATCH --time=10:00:00
#SBATCH --mem=100G
#SBATCH --partition=standard

#1661430
#1661436
#1661454
#1661923
#1813062
#1813100
#1813563
#2066286
echo $(date)
# Load  Julia environment
module load miniforge/24.3.0-py3.11
echo "MODULES LOADED"
python plot_timepoints_critical_radius.py

echo "DONE"