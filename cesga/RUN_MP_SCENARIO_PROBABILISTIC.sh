#!/bin/bash
#SBATCH --job-name=MSE_MP
#SBATCH --output=outputs/mse/logs/MSE_MP_%A_%a.out
#SBATCH --error=outputs/mse/logs/MSE_MP_%A_%a.err

#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G

#SBATCH --array=1-10

set -euo pipefail

module purge
module load cesga/system
module load R/4.4.2

cd "$SLURM_SUBMIT_DIR"

mkdir -p outputs/mse/logs

if [ "$#" -ne 1 ]; then
    echo "Usage: sbatch cesga/RUN_MP_SCENARIO_PROBABILISTIC.sh <scenario_number>"
    exit 1
fi

SCENARIO=$1
BLOCK=$SLURM_ARRAY_TASK_ID

echo "======================================="
echo "PROBABILISTIC FCAP + BESC MSE"
echo "JOB ID: $SLURM_JOB_ID"
echo "TASK ID: $SLURM_ARRAY_TASK_ID"
echo "NODE: $(hostname)"
echo "START: $(date)"
echo "WORKDIR: $(pwd)"
echo "======================================="

echo "SCENARIO: $SCENARIO"
echo "BLOCK: $BLOCK / 10"
echo "======================================="

SCENARIO_TABLE="data/mse/scenarios/refined/candidate_MP_grid_refined_new_RUN.csv"

Rscript scripts/20_RUN_MP_SCENARIO_PROBABILISTIC.R "$SCENARIO" "$BLOCK" "$SCENARIO_TABLE"

echo "======================================="
echo "SCENARIO $SCENARIO - BLOCK $BLOCK COMPLETED"
echo "END: $(date)"
echo "======================================="

