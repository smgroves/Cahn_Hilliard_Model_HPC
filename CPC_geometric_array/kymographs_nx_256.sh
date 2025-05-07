#!/bin/bash
#SBATCH -N 1
#SBATCH -o ../Reports/%A/output.%j.out 
#SBATCH --ntasks-per-node=16
#SBATCH --account=janeslab
#SBATCH --time=10:00:00
#SBATCH --mem=80G
#SBATCH --partition=standard
#SBATCH --array=1-45

export MATLAB_DISABLE_SERVICEHOST=true
module load matlab
indir="/project/g_bme-janeslab/SarahG/julia_out/CPC_geometry/CPC_domain_0_1_e_0.0075"

OPTS=$(sed -n "${SLURM_ARRAY_TASK_ID}"p CPC_geometric_array_eps_0.0075.txt)
dtout=10
dt=0.000001525878906
nx=256
CPC=$(echo "$OPTS" | awk '{print $1}')
cohesin=$(echo "$OPTS" | awk '{print $2}')
time=$(echo "$OPTS" | awk '{print $4}')

echo $CPC
echo $cohesin
echo $time
# Check if the entered number is even
if [ $time == 0.015 ]; then
timesteps=9830
elif [ $time == 0.03 ]; then
timesteps=19661
elif [ $time == 0.045 ]; then
timesteps=29491
fi

echo $timesteps
# if [ $timesteps == 29491 ]; then
name="phi_${nx}_${timesteps}_1.0e-5__CPC_${CPC}_cohesin_${cohesin}_eps_0.0075_domain_0_1"
echo $name 
outdir="$indir/$name"
echo $outdir
matlab -nodisplay -nosplash -r "kymograph_central_droplets_domain('$indir', '$outdir','$name', $dt, $dtout, false, 0, 3.2);quit;"
echo "Done."
# fi

