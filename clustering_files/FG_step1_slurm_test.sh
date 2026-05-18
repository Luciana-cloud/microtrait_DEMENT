#!/bin/bash
#SBATCH --job-name=FG_test_5000
#SBATCH --output=output_%j.txt
#SBATCH --error=error_output_%j.txt
#SBATCH --time=48:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=10
#SBATCH --partition=main
#SBATCH --qos=std
#SBATCH --mail-type=ALL
#SBATCH --mail-user=luciana.chavezrodriguez@wur.nl

# =============================================================================
# SLURM script for Step 1 test run on Anunna

mkdir -p logs
mkdir -p /lustre/nobackup/WUR/ESG/chave013/clustering_YAS/Output_Data/Clustering

echo "================================================="
echo "Job started: $(date)"
echo "Node: $(hostname)"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Memory: 64G"
echo "================================================="

# Load R — check available version first with: module avail R
module purge
module load 2023
ml load R/4.3.2

# Run test script
echo "Starting clustering test..."
Rscript /lustre/nobackup/WUR/ESG/chave013/clustering_YAS/FG_step1_test.R

echo "================================================="
echo "Job finished: $(date)"
echo "================================================="
