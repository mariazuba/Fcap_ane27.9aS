#!/bin/bash
#SBATCH --job-name=MSE_ROB
#SBATCH --output=outputs/mse/logs/MSE_ROB_%A_%a.out
#SBATCH --error=outputs/mse/logs/MSE_ROB_%A_%a.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --array=1-10

set -euo pipefail

# Create outputs/mse/logs BEFORE sbatch, so Slurm can open the logs.
# Usage: sbatch cesga/RUN_MP_ROBUSTNESS.sh <robustness_id> <reference_MP> <prepare|test|full>
# For a single validation block, add --array=1 to sbatch.
if [ "$#" -ne 3 ]; then
    echo "Usage: sbatch cesga/RUN_MP_ROBUSTNESS.sh <robustness_id> <reference_MP> <prepare|test|full>"
    exit 1
fi
ROBUSTNESS_ID=$1
REFERENCE_MP=$2
EXECUTION_MODE=$3
case "$ROBUSTNESS_ID" in REF|R-OM|S-C1|S-C2|I0|I1) ;; *) exit 1 ;; esac
case "$EXECUTION_MODE" in prepare|test|full) ;; *) exit 1 ;; esac
[[ "$REFERENCE_MP" =~ ^[0-9]+$ ]] || exit 1

module purge
module load cesga/system
module load R/4.4.2
cd "$SLURM_SUBMIT_DIR"

printf 'Robustness: %s | Reference MP: %s | Block: %s | Mode: %s\n' "$ROBUSTNESS_ID" "$REFERENCE_MP" "$SLURM_ARRAY_TASK_ID" "$EXECUTION_MODE"
Rscript scripts/27_RUN_MP_ROBUSTNESS.R "$ROBUSTNESS_ID" "$REFERENCE_MP" "$SLURM_ARRAY_TASK_ID" "$EXECUTION_MODE"
