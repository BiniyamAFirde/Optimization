#!/bin/bash -l
#SBATCH --job-name=CPG_k025
#SBATCH --output=Nest_k025_%A_%a.slurmout
#SBATCH --error=Nest_k025_%A_%a.slurmerr
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=0-7
#SBATCH --cpus-per-task=4
#SBATCH --time=04:00:00
#SBATCH --partition=acc

export LANG=${LANG:-C.UTF-8}
export LC_ALL=${LC_ALL:-C.UTF-8}
export PYTHONIOENCODING=utf-8
export PYTHONUNBUFFERED=1

OUTDIR="results_k025_smooth/"
mkdir -p "$OUTDIR"

SEEDS=(12345 22345 32345 42345 52345 62345 72345 82345)
SEED=${SEEDS[${SLURM_ARRAY_TASK_ID:-0}]}

echo "[k025] seed=$SEED  k=0.25 m=0.40 comp=3.0 X0_E=1.10 K_E=15 X0_F=1.42 K_F=30"

srun --cpu-bind=cores \
  python3 -u cpg_final_k025_m040.py \
    --tag        "k025_m040_seed${SEED}" \
    --out        "cpg_run.h5" \
    --outdir     "$OUTDIR" \
    --seed       "$SEED" \
    --sweep-pairs "3.5:0.30" \
    --sweep-run-idx 0 \
    --sweep-dist lognormal_cv \
    --sim-ms     30000 \
    --dt-ms      10 \
    --threads    1 \
    --resolution-ms 0.2 \
    --simulate-chunk-ms 20 \
    --rate-update-ms 20 \
    --weight-sample-ms 1000 \
    --save-weights none \
    --nest-verbosity M_ERROR \
    --delay-model length_velocity \
    --species rat \
    --delay-jitter-ms 0.2 \
    --bs-base-hz 6 \
    --bs-noise-std-hz 0.25 \
    --enforce-tonic-bs \
    --paced-gait \
    --step-period-ms 520 \
    --stance-fraction 0.5 \
    --n-ia-groups 3 \
    --ia-ext-hz 60 80 100 \
    --ia-ext-f-hz 80 \
    --long-run
