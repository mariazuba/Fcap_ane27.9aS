
```bash
#!/bin/bash
#SBATCH --job-name=ane_fcap
#SBATCH --output=outputs/mse/logs/sc_%A_%a.out
#SBATCH --error=outputs/mse/logs/sc_%A_%a.err

#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

#SBATCH --array=1-5

module purge
module load R/4.2.2

cd $SLURM_SUBMIT_DIR

echo "======================================="
echo "JOB ID: $SLURM_JOB_ID"
echo "TASK ID: $SLURM_ARRAY_TASK_ID"
echo "NODE: $(hostname)"
echo "START: $(date)"
echo "======================================="

Rscript scripts/model.R $SLURM_ARRAY_TASK_ID

echo "======================================="
echo "END: $(date)"
echo "======================================="
```
