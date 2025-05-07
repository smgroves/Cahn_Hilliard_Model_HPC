#!/bin/bash
#SBATCH -N 1
#SBATCH -o ../Reports/output.%j.out 
#SBATCH --ntasks-per-node=16
#SBATCH --account=janeslab
#SBATCH --time=10:00:00
#SBATCH --mem=50G
#SBATCH --partition=standard

echo $(date)
module load julia/1.9.2
echo "LOADED JULIA"
# Run Matlab single core program
julia CPC_geometric_run.jl 
echo "DONE"