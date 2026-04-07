#!/bin/bash
#SBATCH -N 1
#SBATCH -o ./Reports/python/%A/output.%j.out 
#SBATCH --ntasks-per-node=16
#SBATCH --account=janeslab
#SBATCH --time=30:00:00
#SBATCH --mem=500G
#SBATCH --partition=standard
#SBATCH --array=1-32

echo $(date)

# Capture initial CPU info
echo "=== CPU SPECIFICATIONS ==="
lscpu | grep -E "Model name|CPU max MHz|CPU min MHz|Socket"
echo ""

# Read the line corresponding to this array task
OPTS=$(sed -n "${SLURM_ARRAY_TASK_ID}p" opts.txt)
echo "Task $SLURM_ARRAY_TASK_ID options: $OPTS"
echo ""

# Parse the options line into variables
# Format: N boundary print method
GridSize=$(echo "$OPTS" | awk '{print $1}')
boundary=$(echo "$OPTS" | awk '{print $2}')
print=$(echo "$OPTS" | awk '{print $3}')
solver=$(echo "$OPTS" | awk '{print $4}')

echo "Parsed parameters:"
echo "  N: $N"
echo "  boundary: $boundary"
echo "  print: $print"
echo "  method: $method"
echo ""

# Load Python environment (adjust module name if needed)
echo $(date)
module load apptainer
module load miniforge
echo "Finished loading modules..."
echo "MODULES LOADED"
echo ""

#64 periodic true NMG


# Run the Python analysis
echo "Running Python analysis..."
python "$./structure_factor_analysis/sf_main.py" ${GridSize} ${boundary} ${print} ${solver} 

echo $(date)
echo "DONE"