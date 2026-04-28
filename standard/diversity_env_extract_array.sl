#!/bin/bash
#SBATCH --job-name=holobiont_sim
#SBATCH --account=uoa04039
#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5G
#SBATCH --array=1-ARRAY_SIZE
#SBATCH --output=holobiont_sim_%A_%a.out
#SBATCH --error=holobiont_sim_%A_%a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=william.pearman@auckland.ac.nz

# Load required modules
module purge
module load R/4.3.1-gimkl-2022a

# Set working directory
cd /nesi/nobackup/uoa04039/Holobiont_PopGen/Simulations_26Nov/simulations_pogen_heatwave/

# Print some job information
echo "Job started at: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Running on node: $SLURM_NODELIST"
echo "CPUs allocated: $SLURM_CPUS_PER_TASK"
echo "Memory allocated: $SLURM_MEM_PER_NODE MB"

# Run R script with array task ID
Rscript ./Analysis_Code/rewritten_analysiscode_array.R $SLURM_ARRAY_TASK_ID

echo "Job finished at: $(date)"
