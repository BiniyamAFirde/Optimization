#!/usr/bin/env bash

set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY="${PY:-python3}"
DRIVER="${DRIVER:-$HERE/cpg_adapter.py}"
GATE_DIAG="${GATE_DIAG:-$HERE/gate_diag.py}"
REPORT="${REPORT:-$HERE/report.py}"
SEEDS="${SEEDS:-12345 22345 32345 42345 52345 62345 72345 82345}"
SIM_MS="${SIM_MS:-30000}"
STOCK_X0F="${STOCK_X0F:-1.42}"
BASELINE="${BASELINE:--0.953}"
FORCE="${FORCE:-0}"

COMMON=(
  --paced-gait --threads 1 --dt-ms 10 --resolution-ms 0.2
  --simulate-chunk-ms 20 --rate-update-ms 20 --long-run
  --nest-verbosity M_ERROR --delay-model length_velocity --species rat
  --delay-jitter-ms 0.2 --bs-base-hz 6 --bs-noise-std-hz 0.25
  --enforce-tonic-bs --stance-fraction 0.5 --n-ia-groups 3
  --ia-ext-hz 60 80 100 --ia-ext-f-hz 80 --static-weight-cv 0.5
  --step-period-ms 520 --inh-comp-f 1.0
  --act-gate-x0-e 1.10 --act-gate-k-e 15
)

usage () {
    echo "usage: run_experiment.sh <grid.tsv> check|passA|plan|passB|report"
    exit 1
}

GRID="${1:-}"
CMD="${2:-}"
[ -n "$GRID" ] && [ -n "$CMD" ] || usage
[ -f "$GRID" ] || { echo "grid not found: $GRID"; exit 1; }
[ -f "$DRIVER" ] || { echo "driver not found: $DRIVER"; exit 1; }

NAME="$(basename "${GRID%.tsv}")"
ROOT="${ROOT:-$HERE/results_$NAME}"

rows () { grep -v '^#' "$GRID" | grep -v '^[[:space:]]*$'; }

run_cell () {
    local label=$1 extra=$2 seed=$3 x0f=$4 dir=$5
    local out="$dir/${label}_seed${seed}.h5"
    mkdir -p "$dir"
    if [ -s "$out" ] && [ "$FORCE" != "1" ]; then
        echo "  skip  ${label}_seed${seed}"
        return 0
    fi
    echo "  run   ${label}_seed${seed}  X0_F=$x0f"
    # shellcheck disable=SC2086
    "$PY" -u "$DRIVER" --out "$out" \
        --experiment "$NAME" --run-label "$label" \
        --run-name "${label}_seed${seed}" \
        --seed "$seed" --sim-ms "$SIM_MS" --act-gate-x0-f "$x0f" \
        "${COMMON[@]}" $extra \
        > "$dir/${label}_seed${seed}.log" 2>&1
    local rc=$?
    [ $rc -eq 0 ] && return 0
    echo "        FAILED rc=$rc"
    grep -E 'error:|unrecognized|Traceback' "$dir/${label}_seed${seed}.log" \
        | tail -3 | sed 's/^/          /'
    return $rc
}

gates () {
    "$PY" - "$ROOT/passA" "$HERE" <<'PY'
import glob, os, re, statistics, sys
sys.path.insert(0, sys.argv[2])
import gate_diag as G
by = {}
for p in sorted(glob.glob(os.path.join(sys.argv[1], "*.h5"))):
    label = re.sub(r"_seed\d+$", "", os.path.splitext(os.path.basename(p))[0])
    try:
        meta, tr, _ = G.load(p, "L")
        r = G.analyse(meta, tr, 10.0)
    except Exception:
        continue
    by.setdefault(label, []).append(r)
for label in sorted(by):
    v = by[label]
    print("%s\t%.3f\t%.1f\t%.3f\t%d" % (
        label,
        statistics.median(x["X0F_rec"] for x in v),
        statistics.median(x["rgf_p50"] for x in v),
        statistics.median(x["corr_RG"] for x in v),
        len(v)))
PY
}

sweep () {
    local dir=$1 mode=$2 table=$3 label extra seed x0f
    while IFS=$'\t' read -r label extra; do
        if [ "$mode" = "own" ]; then
            x0f=$(echo "$table" | awk -F'\t' -v l="$label" '$1==l {print $2}')
            [ -n "$x0f" ] || { echo "  skip  $label (no pass-A gate)"; continue; }
        else
            x0f="$STOCK_X0F"
        fi
        for seed in $SEEDS; do
            run_cell "$label" "$extra" "$seed" "$x0f" "$dir"
        done
    done <<< "$(rows)"
}

case "$CMD" in

check)
    "$PY" -c "import numpy, h5py, matplotlib, nest; print('env ok')" || exit 1
    echo "grid: $GRID"
    rows | awk -F'\t' '{printf "  %-24s %s\n", $1, $2}'
    echo "seeds: $SEEDS"
    echo "runs:  $(( $(rows | wc -l) * $(echo $SEEDS | wc -w) * 2 ))"
    label=$(rows | head -1 | cut -f1)
    extra=$(rows | head -1 | cut -f2)
    SIM_MS=3000 FORCE=1 run_cell "probe" "$extra" 12345 "$STOCK_X0F" \
        "$ROOT/.probe" && echo "probe ok"
    rm -rf "$ROOT/.probe"
    ;;

passA)
    echo "pass A: stock gate $STOCK_X0F"
    sweep "$ROOT/passA" stock ""
    "$PY" "$GATE_DIAG" "$ROOT/passA/*.h5" --csv "$HERE/gate_${NAME}_a.csv"
    ;;

plan)
    t=$(gates)
    [ -n "$t" ] || { echo "no pass-A files in $ROOT/passA"; exit 1; }
    printf '%-26s %10s %10s %9s %6s\n' condition "RG-F p50" "gate B" corr_RG seeds
    echo "$t" | while IFS=$'\t' read -r l x0 p50 crg n; do
        printf '%-26s %10s %10s %9s %6s\n' "$l" "$p50" "$x0" "$crg" "$n"
    done
    ;;

passB)
    t=$(gates)
    [ -n "$t" ] || { echo "run passA first"; exit 1; }
    echo "pass B: each condition at its own gate"
    sweep "$ROOT/passB" own "$t"
    "$PY" "$GATE_DIAG" "$ROOT/passB/*.h5" --csv "$HERE/gate_${NAME}_b.csv"
    ;;

report)
    "$PY" "$REPORT" "$HERE/gate_${NAME}_b.csv" --baseline "$BASELINE"
    ;;

*)
    usage
    ;;
esac
