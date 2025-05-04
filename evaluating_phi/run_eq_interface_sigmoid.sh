#!/bin/bash
#SBATCH -N 1
#SBATCH -o ./Reports/%A/output.%j.out 
#SBATCH --ntasks-per-node=16
#SBATCH --account=janeslab
#SBATCH --time=10:00:00
#SBATCH --mem=100G
#SBATCH --partition=standard
module load matlab

matlab -nodisplay -nosplash -r "flat_interface_eq_interface_sigmoid();quit;"

