#!/bin/bash -e

#SBATCH --job-name=holobiont_sim
#SBATCH --account=uoa04039
#SBATCH --time=4:00:00              # Adjust based on your simulation time
#SBATCH --mem=3G                     # Memory per task
#SBATCH --cpus-per-task=1            # 1 CPU per simulation
#SBATCH --array=1-720  # Replace ARRAY_MAX with actual number; %200 limits concurrent jobs
#SBATCH --output=logs/holobiont_%A_%a.out
#SBATCH --error=logs/holobiont_%A_%a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=william.pearman@auckland.ac.nz

# =====================================================================
# HOLOBIONT SIMULATION - SLURM ARRAY JOB
# Each array task runs one simulation (combination + replicate)
# =====================================================================

# Print diagnostic information
echo "=========================================="
echo "SLURM JOB INFORMATION"
echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Node: $SLURM_NODELIST"
echo "Job name: $SLURM_JOB_NAME"
echo "Start time: $(date)"
echo "=========================================="
echo ""

# Load required modules
module purge
module load R/4.3.1-gimkl-2022a  # Adjust to your R version

# Set working directory
cd /nesi/nobackup/uoa04039/Holobiont_PopGen/Simulations_26Nov/simulations_pogen_heatwave/Rerun_sims_array_28Nov/Importance_Mutation

# Create logs directory if it doesn't exist
mkdir -p logs

# Run the simulation
echo "Starting simulation for array task $SLURM_ARRAY_TASK_ID"
echo ""

Rscript holobiont_simulation_array.R $SLURM_ARRAY_TASK_ID

# Check exit status
exit_status=$?
if [ $exit_status -eq 0 ]; then
    echo ""
    echo "=========================================="
    echo "Task $SLURM_ARRAY_TASK_ID completed successfully"
    echo "End time: $(date)"
    echo "=========================================="
else
    echo ""
    echo "=========================================="
    echo "Task $SLURM_ARRAY_TASK_ID FAILED with exit status $exit_status"
    echo "End time: $(date)"
    echo "=========================================="
fi

exit $exit_status
