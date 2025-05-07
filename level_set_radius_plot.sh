#!/bin/bash
#SBATCH -N 1
#SBATCH -o ./Reports/output.%j.out 
#SBATCH --ntasks-per-node=16
#SBATCH --account=janeslab
#SBATCH --time=10:00:00
#SBATCH --mem=50G
#SBATCH --partition=standard

echo $(date)
module load matlab
echo "LOADED MATLAB"
# Run Matlab single core program
matlab -nodisplay -r "level_set_radius(); exit;"
echo "DONE"
