#!/bin/bash
#SBATCH -N 1
#SBATCH -o ../Reports/output.%j.out 
#SBATCH --ntasks-per-node=16
#SBATCH --account=janeslab
#SBATCH --time=26:00:00
#SBATCH --mem=400G
#SBATCH --partition=standard

module load matlab

indir="/project/g_bme-janeslab/SarahG/julia_out/CPC_geometry/CPC_domain_0_2_e_0.0067_t_0.4"
name="phi_512_262144_1.0e-5__CPC_0.12_cohesin_0.14_eps_0.0067_domain_0_2"
echo $name
matlab -nodisplay -nosplash -r "size_phi('$indir', '$name');quit;"

name="phi_512_262144_1.0e-5__CPC_0.35_cohesin_0.14_eps_0.0067_domain_0_2"
echo $name
matlab -nodisplay -nosplash -r "size_phi('$indir', '$name');quit;"

name="phi_512_262144_1.0e-5__CPC_0.35_cohesin_0.12_eps_0.0067_domain_0_2"
echo $name
matlab -nodisplay -nosplash -r "size_phi('$indir', '$name');quit;"

#4133015