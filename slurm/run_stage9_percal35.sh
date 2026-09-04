#!/bin/bash -l
#SBATCH --job-name=CPG_s9
#SBATCH --output=Nest_s9_%A_%a.slurmout
#SBATCH --error=Nest_s9_%A_%a.slurmerr
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=0-7%5
#SBATCH --cpus-per-task=4
#SBATCH --time=04:00:00
#SBATCH --partition=acc

export LANG=${LANG:-C.UTF-8}
export LC_ALL=${LC_ALL:-C.UTF-8}
export PYTHONIOENCODING=utf-8
export PYTHONUNBUFFERED=1

GATES="gates_k020_m035.tsv"
if [ ! -s "$GATES" ]; then
    echo "ERROR: $GATES missing or empty."
    echo "Run:  python3 make_gate_table.py \"results_stage8/*.h5\" > gates.tsv"
    exit 1
fi

OUTDIR="results_stage9/"
mkdir -p "$OUTDIR"

T=${SLURM_ARRAY_TASK_ID:-0}
K=0.20; M=0.35; COMP=20.0

# Row T+1 of the table 
LINE=$(sed -n "$((T+1))p" "$GATES")
if [ -z "$LINE" ]; then echo "No row $((T+1)) in $GATES"; exit 1; fi
SEED=$(echo "$LINE" | cut -f1)
X0E=$(echo "$LINE"  | cut -f2)
KE=$(echo "$LINE"   | cut -f3)
X0F=$(echo "$LINE"  | cut -f4)
KF=$(echo "$LINE"   | cut -f5)

echo "[PERCAL] task=$T seed=$SEED threads=1 k=$K m=$M comp=$COMP X0_E=$X0E K_E=$KE X0_F=$X0F K_F=$KF"

srun --cpu-bind=cores \
  python3 -u cpg_variantAF_seeded.py \
    --tag        "percal_s8_seed${SEED}" \
    --out        "cpg_run.h5" \
    --outdir     "$OUTDIR" \
    --seed       "$SEED" \
    --k-conn     "$K" \
    --m-neurons  "$M" \
    --inh-comp   "$COMP" \
    --inh-comp-f 1.0 \
    --act-gate-x0-e "$X0E" \
    --act-gate-k-e  "$KE" \
    --act-gate-x0-f "$X0F" \
    --act-gate-k-f  "$KF" \
    --sweep-pairs "3.5:0.30" \
    --sweep-run-idx 0 \
    --sweep-dist lognormal_cv \
    --sim-ms     30000 \
    --dt-ms      10 \
    --threads    1 \
    --resolution-ms 0.2 \
    --simulate-chunk-ms 100 \
    --rate-update-ms 100 \
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
