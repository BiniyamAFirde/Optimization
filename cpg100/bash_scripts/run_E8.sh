#!/bin/bash -l
#SBATCH --job-name=cpg_E8
#SBATCH --output=slurm_E8_%A_%a.out
#SBATCH --error=slurm_E8_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=01:30:00
#SBATCH --mem=2G
#
# run_E8.sh -- roadmap node E8: final validation of the 100-neuron recipes
# =======================================================================
#
#   bash run_E8.sh check             files + python packages + NEST present?  (run first)
#   bash run_E8.sh plan              cells, task count (no NEST)
#   bash run_E8.sh counts            build each network once, print ACTUAL neuron/synapse counts
#   bash run_E8.sh smoke             5 s, n100 STDP recipe at 350 ms, pass A then B
#   bash run_E8.sh --local [A|B|all] everything serially (Mac)   ~2.5-3 h (baseline is 2/3 of it)
#   bash run_E8.sh submit            Slurm: pass-A array, pass-B array (afterany)
#   bash run_E8.sh summarize         gait-by-gait verdict, both legs, L-R phase (pass B)
#
#   Staging on the Mac (100-neuron cells first, ~25 min; then the baseline):
#     CELLS_ONLY="$(bash run_E8.sh list n100)" bash run_E8.sh --local all
#     bash run_E8.sh --local all
#
# ── WHY ──────────────────────────────────────────────────────────────────
#   E6 (Mac, 60 s, 8 seeds) confirmed the 100-neuron recipes. Two finalists:
#     A  n100_kmin2c, STDP off     -0.949 +- 0.027 (baseline off -0.953), trough 1.32,
#                                  2,353 synapses onto neurons
#     B  n100_kmin3c_inex15, STDP 1e-6
#                                  -0.974 +- 0.015 (baseline on -0.973), trough 0.91,
#                                  weights near-settled (2.8 %/min; 1e-5 drifts 22 %/min)
#   Every result so far is leg L at 520 ms. E8 checks what the roadmap requires
#   before calling it done: BOTH legs, left-right phase, and gaits 350 / 520 / 780 ms,
#   each against the k=1 baseline run at the SAME gait.
#
# ── GRID: 9 cells x 8 seeds x 2 passes = 144 runs, 60 s each ─────────────
#   base_k100_g{350,520,780}_stdpoff           k=1 baseline, 1200 neurons (reference per gait)
#   n100_kmin2c_g{350,520,780}_stdpoff         finalist A
#   n100_kmin3c_inex15_lam6_g{..}_stdpon       finalist B
#   Per-leg two-pass gate, chunk 100 ms, 1 thread.
#
# ── PASS (summarize --by-gait, summarize_small v1.6) ──────────────────────
#   A1 |corr_F - baseline(same gait)| <= 0.02   (STDP-on: baseline - 0.020, E0 shift)
#   A2 troughF <= max(4, 1.25 x baseline trough at that gait)   A3 SD <= 0.03
#   A4 readout not under-reporting   A5 stable over 60 s
#   A6 RIGHT leg within 0.02 of the baseline's right leg   A7 L-R phase 0.5 +- 0.05, R >= 0.9
#
# ── WHAT DECIDES THE NEXT NODE ───────────────────────────────────────────
#   a finalist passes at all three gaits -> DONE: write up (FINDINGS + final figure)
#   fails at one gait                    -> E3 at that gait (inhibition tuning)

set -uo pipefail

HERE="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
PY="${PY:-python3}"
DRIVER="${DRIVER:-$HERE/cpg_e8.py}"
SUMMARY="${SUMMARY:-$HERE/summarize_small.py}"
ROOT="${ROOT:-$HERE/results_E8}"
SEEDS="${SEEDS:-12345 22345 32345 42345 52345 62345 72345 82345}"
SIM_S="${SIM_S:-60}"
LAMBDA="${LAMBDA:-1e-5}"
PARTITION="${PARTITION:-acc}"
TIME="${TIME:-01:30:00}"
FORCE="${FORCE:-0}"
STRICT_NEST="${STRICT_NEST:-0}"
X0F_STOCK=1.42
# k=1 baseline trace for the top panel of the figure (from E0), if present
BASELINE_H5="${BASELINE_H5:-$(ls "$HERE"/results_E0/passB/base_k100_stdpoff_seed*.h5 2>/dev/null | head -1)}"

# name | preset | stdp | rule | preserve(0/1) | extra flags (space-separated, may be empty)
# 100-neuron cells first so a partial run is already useful; the baseline takes ~3 min/run.
CELLS=(
  "n100_kmin2c_g350_stdpoff|n100|off|indegree|0|--indeg-min 2 --indeg-min-conserve --step-period-ms 350"
  "n100_kmin2c_g520_stdpoff|n100|off|indegree|0|--indeg-min 2 --indeg-min-conserve --step-period-ms 520"
  "n100_kmin2c_g780_stdpoff|n100|off|indegree|0|--indeg-min 2 --indeg-min-conserve --step-period-ms 780"
  "n100_kmin3c_inex15_lam6_g350_stdpon|n100|on|indegree|0|--indeg-min 3 --indeg-min-conserve --inh-comp 5.625 --stdp-lambda 1e-6 --step-period-ms 350"
  "n100_kmin3c_inex15_lam6_g520_stdpon|n100|on|indegree|0|--indeg-min 3 --indeg-min-conserve --inh-comp 5.625 --stdp-lambda 1e-6 --step-period-ms 520"
  "n100_kmin3c_inex15_lam6_g780_stdpon|n100|on|indegree|0|--indeg-min 3 --indeg-min-conserve --inh-comp 5.625 --stdp-lambda 1e-6 --step-period-ms 780"
  "base_k100_g350_stdpoff|base_k100|off|bernoulli|0|--step-period-ms 350"
  "base_k100_g520_stdpoff|base_k100|off|bernoulli|0|--step-period-ms 520"
  "base_k100_g780_stdpoff|base_k100|off|bernoulli|0|--step-period-ms 780"
)
say () { printf '\n\033[1m== %s\033[0m\n' "$*"; }
die () { printf '\n[FATAL] %s\n' "$*" >&2; exit 1; }

if [ -n "${CELLS_ONLY:-}" ]; then
  _keep=()
  for c in "${CELLS[@]}"; do
    for w in $CELLS_ONLY; do [ "${c%%|*}" = "$w" ] && _keep+=("$c"); done
  done
  [ ${#_keep[@]} -gt 0 ] || die "CELLS_ONLY='$CELLS_ONLY' matches no cell (see: bash run_E8.sh plan)"
  CELLS=("${_keep[@]}")
fi
read -r -a SEED_ARR <<< "$SEEDS"
NCELL=${#CELLS[@]}
NSEED=${#SEED_ARR[@]}
NTASK=$((NCELL * NSEED))

if [ "${1:-}" != "check" ]; then
    for f in cpg_small.py cpg_e0.py cpg_e1.py cpg_e2.py cpg_e8.py; do
        [ -f "$HERE/$f" ] || die "$f must be in $HERE -- run: bash run_E8.sh check"
    done
fi
cell_field () { echo "$1" | cut -d'|' -f"$2"; }

# prescribe <cell> -> "X0F_L X0F_R" for pass B (per-leg median over pass-A seeds)
prescribe () {
    "$PY" - "$ROOT/passA" "$1" <<'PY'
import glob, json, os, statistics, sys
d, cell = sys.argv[1], sys.argv[2]
L, R = [], []
for p in sorted(glob.glob(os.path.join(d, f"{cell}_seed*.json"))):
    if p.endswith(".counts.json"):
        continue
    try:
        r = json.load(open(p))
        L.append(r["legs"]["leg_L"]["X0F_rec"]); R.append(r["legs"]["leg_R"]["X0F_rec"])
    except Exception as e:
        print(f"[prescribe] skip {os.path.basename(p)}: {e!r}", file=sys.stderr)
if not L:
    sys.exit(f"[prescribe] no pass-A results for {cell} in {d}")
xl, xr = statistics.median(L), statistics.median(R)
print(f"[prescribe] {cell}: X0_F L = {xl:.3f}  R = {xr:.3f}  (median of {len(L)} seeds)",
      file=sys.stderr)
print(f"{xl:.3f} {xr:.3f}")
PY
}

# run_one <A|B> <cell-spec> <seed> <x0f_L> <x0f_R>
run_one () {
    local pass=$1 spec=$2 seed=$3 xl=$4 xr=$5
    local cell preset stdp rule pres
    cell=$(cell_field "$spec" 1); preset=$(cell_field "$spec" 2); stdp=$(cell_field "$spec" 3)
    rule=$(cell_field "$spec" 4); pres=$(cell_field "$spec" 5)
    local outdir="$ROOT/pass$pass" tag="${cell}_seed${seed}"
    mkdir -p "$outdir"
    if [ -s "$outdir/$tag.json" ] && [ "$FORCE" != "1" ]; then
        echo "  [skip] pass $pass $tag"; return 0
    fi
    local extra=() more
    [ "$pres" = "1" ] && extra+=(--preserve-input)
    [ "$STRICT_NEST" = "1" ] && extra+=(--strict-nest)
    more=$(cell_field "$spec" 6)
    if [ -n "$more" ]; then read -r -a _m <<< "$more"; extra+=("${_m[@]}"); fi
    echo "  [run ] pass $pass $tag  stdp=$stdp ${more:-(no extra flags)} X0_F L/R=$xl/$xr  (${SIM_S}s)"
    "$PY" -u "$DRIVER" --preset "$preset" --out "$outdir" --tag "$tag" --cell "$cell" \
        --pass-label "$pass" --stdp "$stdp" --stdp-lambda "$LAMBDA" --conn-rule "$rule" \
        --sim-s "$SIM_S" --seed "$seed" --x0f "$xl" --x0f-r "$xr" --x0e 1.10 --threads 1 \
        ${extra[@]+"${extra[@]}"} > "$outdir/$tag.log" 2>&1
    local rc=$?
    if [ $rc -ne 0 ]; then
        echo "        FAILED rc=$rc -- tail of $outdir/$tag.log:"
        tail -5 "$outdir/$tag.log" | sed 's/^/        | /'
        return $rc
    fi
    grep -E '^\[(done|E0)\]' "$outdir/$tag.log" | cut -c1-170 | sed 's/^/        /'
}

task () {
    local pass=$1 idx=$2
    [ "$idx" -lt "$NTASK" ] || die "task index $idx >= $NTASK"
    local spec="${CELLS[$((idx / NSEED))]}" seed="${SEED_ARR[$((idx % NSEED))]}"
    local cell; cell=$(cell_field "$spec" 1)
    local xl=$X0F_STOCK xr=$X0F_STOCK x
    if [ "$pass" = "B" ]; then
        x=$(prescribe "$cell") || die "no pass-A prescription for $cell"
        xl=${x% *}; xr=${x#* }
        mkdir -p "$ROOT/passB"; echo "$x" > "$ROOT/passB/$cell.x0f"
    fi
    run_one "$pass" "$spec" "$seed" "$xl" "$xr"
}

cmd_check () {
    say "Pre-flight: $(uname -s) $(uname -m), bash $BASH_VERSION, python = $(command -v "$PY")"
    local bad=0 f
    for f in cpg_small.py cpg_e0.py cpg_e1.py cpg_e2.py cpg_e8.py run_E8.sh summarize_small.py; do
        if [ -f "$HERE/$f" ]; then printf '  ok       %s\n' "$f"
        else printf '  MISSING  %s\n' "$f"; bad=1; fi
    done
    if [ -f "$HERE/cpg_small.py" ] && ! grep -q 'cpg_small 1.2' "$HERE/cpg_small.py"; then
        echo "  OLD      cpg_small.py is not v1.2 -- replace it (adds --indeg-min)"; bad=1
    fi
    if [ -n "$BASELINE_H5" ]; then echo "  ok       baseline trace for figures: ${BASELINE_H5#$HERE/}"
    else echo "  note     no results_E0/passB/base_k100_stdpoff_*.h5 here: figure will have no baseline panel"; fi
    echo "  conda env: ${CONDA_DEFAULT_ENV:-none}"
    PYNEST_QUIET=1 "$PY" - <<'PY' || bad=1
import sys
ok = True
for mod in ("numpy", "h5py", "matplotlib"):
    try:
        m = __import__(mod); print(f"  ok       {mod} {m.__version__}")
    except Exception as e:
        ok = False; print(f"  MISSING  {mod}  ({e.__class__.__name__})")
try:
    import nest
    nest.set_verbosity("M_ERROR")
    print(f"  ok       nest {nest.__version__}")
except Exception as e:
    ok = False; print(f"  MISSING  nest  ({e.__class__.__name__}: {e})")
sys.exit(0 if ok else 1)
PY
    if [ "$bad" = 0 ]; then echo; echo "  All good. Next: bash run_E8.sh smoke"
    else echo; echo "  Not ready: conda activate nest, and put the files listed above in $HERE"; return 1; fi
}

cmd_plan () {
    say "E8 grid: $NCELL cells x $NSEED seeds = $NTASK tasks per pass (x2 passes), ${SIM_S}s each"
    local c; for c in "${CELLS[@]}"; do
        printf '  %-28s preset=%-9s stdp=%-4s %s %s\n' "$(cell_field "$c" 1)" \
            "$(cell_field "$c" 2)" "$(cell_field "$c" 3)" "$(cell_field "$c" 4)" "$(cell_field "$c" 6)"
    done
    echo "  lambda=$LAMBDA  seeds: $SEEDS   results -> $ROOT/passA, $ROOT/passB"
}

cmd_counts () {
    say "Actual counts NEST builds (build only, seed ${SEED_ARR[0]})"
    mkdir -p "$ROOT/counts"
    local c; for c in "${CELLS[@]}"; do
        [ "$(cell_field "$c" 3)" = "off" ] || continue
        local n; n=$(cell_field "$c" 1)
        local extra=() more; [ "$(cell_field "$c" 5)" = "1" ] && extra+=(--preserve-input)
        more=$(cell_field "$c" 6)
        if [ -n "$more" ]; then read -r -a _m <<< "$more"; extra+=("${_m[@]}"); fi
        "$PY" "$DRIVER" --preset "$(cell_field "$c" 2)" --conn-rule "$(cell_field "$c" 4)" \
            ${extra[@]+"${extra[@]}"} --build-only --out "$ROOT/counts" --tag "$n" \
            --seed "${SEED_ARR[0]}" --threads 1 > "$ROOT/counts/$n.build.log" 2>&1 \
            || { echo "  $n: BUILD FAILED"; tail -4 "$ROOT/counts/$n.build.log" | sed 's/^/     | /'; continue; }
        grep -E '^\[counts\] (Izh|syn)' "$ROOT/counts/$n.build.log" | sed "s/^\[counts\]/  $n:/"
    done
}

cmd_local () {
    local which=${1:-all} i
    if [ -z "${E8_AWAKE:-}" ] && command -v caffeinate >/dev/null 2>&1; then
        echo "  (macOS: re-running under caffeinate -i so the Mac does not sleep mid-run)"
        E8_AWAKE=1 exec caffeinate -i bash "$HERE/run_E8.sh" --local "$which"
    fi
    echo "  Serial run, 60 s each: ~15 s wall per 100-neuron run, ~3 min per k=1 baseline run."
    echo "  Finished runs are skipped on restart: Ctrl-C is safe, re-run the same command to resume."
    if [ "$which" = "A" ] || [ "$which" = "all" ]; then
        say "Pass A (stock gate X0_F=$X0F_STOCK, both legs) -- $NTASK runs"
        for ((i = 0; i < NTASK; i++)); do task A "$i"; done
    fi
    if [ "$which" = "B" ] || [ "$which" = "all" ]; then
        say "Pass B (per-cell, per-leg prescribed X0_F) -- $NTASK runs"
        for ((i = 0; i < NTASK; i++)); do task B "$i"; done
    fi
    [ "$which" = "all" ] && cmd_summarize
}

cmd_submit () {
    command -v sbatch >/dev/null || die "sbatch not found -- use: bash run_E8.sh --local"
    mkdir -p "$ROOT"
    local exp="ALL,SEEDS=$SEEDS,SIM_S=$SIM_S,ROOT=$ROOT,LAMBDA=$LAMBDA,PY=$PY"
    exp="$exp,STRICT_NEST=$STRICT_NEST,CELLS_ONLY=${CELLS_ONLY:-}"
    local ja jb
    ja=$(sbatch --parsable -p "$PARTITION" --time="$TIME" --array=0-$((NTASK - 1)) \
         --export="$exp" "$HERE/run_E8.sh" task A) || die "pass-A submit failed"
    jb=$(sbatch --parsable -p "$PARTITION" --time="$TIME" --array=0-$((NTASK - 1)) \
         --export="$exp" --dependency=afterany:"$ja" "$HERE/run_E8.sh" task B) \
         || die "pass-B submit failed"
    echo "  pass A: job $ja   pass B: job $jb   ->  when B is done: bash run_E8.sh summarize"
}

cmd_smoke () {
    local spec="n100_kmin3c_inex15_lam6_g350_stdpon|n100|on|indegree|0|--indeg-min 3 --indeg-min-conserve --inh-comp 5.625 --stdp-lambda 1e-6 --step-period-ms 350" x
    say "Smoke test: n100_kmin3c_inex15_lam6_g350_stdpon, 5 s, seed ${SEED_ARR[0]}, pass A then B"
    ROOT="$ROOT/smoke" SIM_S=5 FORCE=1
    mkdir -p "$ROOT"; run_one A "$spec" "${SEED_ARR[0]}" "$X0F_STOCK" "$X0F_STOCK" || die "pass A failed"
    x=$(prescribe n100_kmin3c_inex15_lam6_g350_stdpon) || die "prescription failed"
    run_one B "$spec" "${SEED_ARR[0]}" "${x% *}" "${x#* }" || die "pass B failed"
    "$PY" "$SUMMARY" --root "$ROOT" --ref-prefix ctrl_n480 --prefix "$ROOT/E8_summary" >/dev/null \
        && echo "  summarizer ok" || echo "  (summarizer FAILED)"
    echo "  ok -- outputs in $ROOT   (5 s, 1 seed: NOT evidence)"
}

cmd_summarize () {
    [ -f "$SUMMARY" ] || die "summarize_small.py not found: $SUMMARY"
    "$PY" "$SUMMARY" --root "$ROOT" --ref-prefix base_k100 --prefix "$ROOT/E8_summary" \
        --t0 40 --by-gait ${BASELINE_H5:+--baseline "$BASELINE_H5"}
}

case "${1:-}" in
    check)           cmd_check ;;
    plan)            cmd_plan ;;
    counts)          cmd_counts ;;
    smoke)           cmd_smoke ;;
    submit)          cmd_submit ;;
    --local|local)   cmd_local "${2:-all}" ;;
    task)            task "${2:?A or B}" "${SLURM_ARRAY_TASK_ID:?not in a Slurm array}" ;;
    summarize)       cmd_summarize ;;
    list)            for c in "${CELLS[@]}"; do n=${c%%|*}; case "$n" in *"${2:-}"*) printf '%s ' "$n";; esac; done; echo ;;
    *) sed -n '13,21p' "${BASH_SOURCE[0]}"; exit 1 ;;
esac
