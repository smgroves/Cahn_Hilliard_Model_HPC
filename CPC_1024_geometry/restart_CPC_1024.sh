#!/bin/bash
#SBATCH -N 1
#SBATCH -o ./Reports/%A/output.%J.out
#SBATCH --ntasks-per-node=16
#SBATCH --account=janeslab
#SBATCH --time=10:00:00
#SBATCH --mem=100G
#SBATCH --partition=standard
echo $(date)

#append the results from the last failed run to the main file.
minimum=19660
search_dir=/scratch/xpz5km/Cahn_Hilliard_Model/julia_out/CPC_geometry_1024
for entry in "$search_dir"/*
do 
    number=$(echo "$entry" | grep -oP '(?<=phi_1024_)\d+')
    echo $number
    if [ "$number" -lt $minimum ]; then 
        minimum=$number
    fi
done

oldfile="/scratch/xpz5km/Cahn_Hilliard_Model/julia_out/CPC_geometry_1024/phi_1024_19660_0.0001__CPC_40_cohesin_16_eps_0.007504684956431058.txt"
newfile="/scratch/xpz5km/Cahn_Hilliard_Model/julia_out/CPC_geometry_1024/phi_1024_${minimum}_0.0001__CPC_40_cohesin_16_eps_0.007504684956431058.txt"
cat $newfile >> $oldfile; rm $newfile

module load anaconda
module load julia/1.9.2
echo "MODULES LOADED"

python last_n_lines.py
julia restart_CPC_1024.jl
echo "DONE"




