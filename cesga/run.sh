#!/bin/bash
#SBATCH --job-name=ane_fcap_30y_1000i
#SBATCH --output=outputs/mse/logs/sc_%A_%a.out
#SBATCH --error=outputs/mse/logs/sc_%A_%a.err

#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G

#SBATCH --array=1-5

module purge
module load R/4.2.2

cd $SLURM_SUBMIT_DIR

mkdir -p outputs/mse/logs
mkdir -p outputs/mse/res
mkdir -p outputs/mse/summary

echo "======================================="
echo "JOB ID: $SLURM_JOB_ID"
echo "TASK ID: $SLURM_ARRAY_TASK_ID"
echo "NODE: $(hostname)"
echo "START: $(date)"
echo "WORKDIR: $(pwd)"
echo "======================================="

Rscript scripts/model.R $SLURM_ARRAY_TASK_ID

echo "======================================="
echo "END: $(date)"
echo "======================================="
