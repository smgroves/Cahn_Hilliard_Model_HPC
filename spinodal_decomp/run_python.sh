#!/bin/bash
#SBATCH -N 1
#SBATCH -o ./Reports/python/%A/output.%j.out 
#SBATCH --ntasks-per-node=16
#SBATCH --account=janeslab
#SBATCH --time=25:00:00
#SBATCH --mem=100G
#SBATCH --partition=standard
#SBATCH --array=1-4

echo $(date)
module load apptainer

#install conda environment if you haven't already
# conda create -n vcell python=3.9
# conda activate vcell
# pip install pandas h5py

echo "Loading modules..."
module load miniforge
echo "Finished loading modules..."

OPTS=$(sed -n "${SLURM_ARRAY_TASK_ID}"p opts_redo_python_v3.txt)
echo $OPTS

outdir="/project/g_bme-janeslab/SarahG/spinodal_decomp_04_2025/out_python"
mkdir -p $outdir

GridSize=$(echo "$OPTS" | awk '{print $1}')
boundary=$(echo "$OPTS" | awk '{print $2}')
print=$(echo "$OPTS" | awk '{print $3}')
solver=$(echo "$OPTS" | awk '{print $4}')

SLURM_ID=${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
python ./CahnHilliard_Python_solvers/run_spinodal_decomp.py ${GridSize} ${boundary} ${print} ${solver} ${SLURM_ID} ""


echo "DONE"
