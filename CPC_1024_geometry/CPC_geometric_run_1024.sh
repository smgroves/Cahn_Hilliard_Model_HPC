#!/bin/bash
#SBATCH -N 1
#SBATCH -o ./Reports/output.%J.out
#SBATCH --ntasks-per-node=16
#SBATCH --account=janeslab
#SBATCH --time=12:00:00
#SBATCH --mem=100G
#SBATCH --partition=standard
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-type=fail         # send email if job fails
#SBATCH --mail-user=xpz5km@virginia.edu

echo $(date)
module load julia/1.9.2
echo "LOADED JULIA"
# Run Matlab single core program
julia run_solver_manuscript_1024.jl 
echo "DONE"