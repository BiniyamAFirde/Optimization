#!/bin/bash -l
#SBATCH --job-name=cpg_E3
#SBATCH --output=slurm_E3_%A_%a.out
#SBATCH --error=slurm_E3_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH --mem=2G
#
# run_E3.sh -- roadmap node E3: connection floor + inhibition tuning at 100 neurons
# ================================================================================
#
#   bash run_E3.sh check             files + python packages + NEST present?  (run first)
#   bash run_E3.sh plan              cells, task count (no NEST)
#   bash run_E3.sh counts            build each network once, print ACTUAL neuron/synapse counts
#   bash run_E3.sh smoke             5 s, n100_kmin3c_stdpoff, pass A then B
#   bash run_E3.sh --local [A|B|all] everything serially (Mac)        ~20-30 min in total
#   bash run_E3.sh submit            Slurm: pass-A array, pass-B array (afterany)
#   bash run_E3.sh summarize         table + symptoms + traces (pass B)
#
# ── WHY ──────────────────────────────────────────────────────────────────
#   E2 (Mac, 8 seeds): fixed in-degree holds to 140 neurons (-0.938 +- 0.019);
#   at 100 it gives -0.914 +- 0.046 with corr_RG -0.911 (the rhythm itself, not
#   the readout) and a deep trough (2.32). At 100 neurons EVERY neuron->neuron
#   pathway has exactly 1 input per target, for any k <= 0.6 -- each neuron's
#   rhythm hangs on one randomly chosen partner. So E3 first raises that floor
#   (cpg_small v1.2 --indeg-min), then tunes inhibition on the best floor.
#
# ── GRID: 10 cells x 8 seeds x 2 passes = 160 runs, 30 s each ────────────
#   cell                        flags on top of --preset n100 --conn-rule indegree
#   ctrl_n480_stdpoff           (480-neuron control, bernoulli)
#   n140_indeg_stdpoff          (E2's smallest passing rung, anchor)
#   n100_indeg_stdpoff          none                         E2 reference, K = 1
#   n100_kmin2c_stdpoff         --indeg-min 2 --indeg-min-conserve
#   n100_kmin3c_stdpoff         --indeg-min 3 --indeg-min-conserve
#   n100_kmin3_stdpoff          --indeg-min 3        (no conserve: 3x coupling)
#   n100_kmin3c_infx2_stdpoff   kmin3c + --inh-comp-f 2.0   (F->InF->E x2: silence E harder)
#   n100_kmin3c_inex05_stdpoff  kmin3c + --inh-comp 1.875   (E->InE->F x0.5)
#   n100_kmin3c_inex15_stdpoff  kmin3c + --inh-comp 5.625   (E->InE->F x1.5)
#   n100_kmin3c_stdpon          kmin3c + STDP (lambda 1e-5)
#   n100 = 100 Izhikevich neurons; k = 0.20; chunk 100 ms; per-leg two-pass gate.
#
# ── WHAT DECIDES THE NEXT NODE ───────────────────────────────────────────
#   a n100 cell passes                   -> E6 (STDP arm, 60 s) with that recipe
#   floor helps SD but corr_RG stays low -> E4 (drive) / E5 (intrinsic RG bursting)
#   nothing helps                        -> E4

set -uo pipefail

HERE="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
PY="${PY:-python3}"
DRIVER="${DRIVER:-$HERE/cpg_e3.py}"
SUMMARY="${SUMMARY:-$HERE/summarize_small.py}"
ROOT="${ROOT:-$HERE/results_E3}"
SEEDS="${SEEDS:-12345 22345 32345 42345 52345 62345 72345 82345}"
SIM_S="${SIM_S:-30}"
LAMBDA="${LAMBDA:-1e-5}"
PARTITION="${PARTITION:-acc}"
TIME="${TIME:-00:30:00}"
FORCE="${FORCE:-0}"
STRICT_NEST="${STRICT_NEST:-0}"
X0F_STOCK=1.42
# k=1 baseline trace for the top panel of the figure (from E0), if present
BASELINE_H5="${BASELINE_H5:-$(ls "$HERE"/results_E0/passB/base_k100_stdpoff_seed*.h5 2>/dev/null | head -1)}"

# name | preset | stdp | rule | preserve(0/1) | extra flags (space-separated, may be empty)
CELLS=(
  "ctrl_n480_stdpoff|ctrl_n480|off|bernoulli|0|"
  "n140_indeg_stdpoff|n140|off|indegree|0|"
  "n100_indeg_stdpoff|n100|off|indegree|0|"
  "n100_kmin2c_stdpoff|n100|off|indegree|0|--indeg-min 2 --indeg-min-conserve"
  "n100_kmin3c_stdpoff|n100|off|indegree|0|--indeg-min 3 --indeg-min-conserve"
  "n100_kmin3_stdpoff|n100|off|indegree|0|--indeg-min 3"
  "n100_kmin3c_infx2_stdpoff|n100|off|indegree|0|--indeg-min 3 --indeg-min-conserve --inh-comp-f 2.0"
  "n100_kmin3c_inex05_stdpoff|n100|off|indegree|0|--indeg-min 3 --indeg-min-conserve --inh-comp 1.875"
  "n100_kmin3c_inex15_stdpoff|n100|off|indegree|0|--indeg-min 3 --indeg-min-conserve --inh-comp 5.625"
  "n100_kmin3c_stdpon|n100|on|indegree|0|--indeg-min 3 --indeg-min-conserve"
)
say () { printf '\n\033[1m== %s\033[0m\n' "$*"; }
die () { printf '\n[FATAL] %s\n' "$*" >&2; exit 1; }

if [ -n "${CELLS_ONLY:-}" ]; then
  _keep=()
  for c in "${CELLS[@]}"; do
    for w in $CELLS_ONLY; do [ "${c%%|*}" = "$w" ] && _keep+=("$c"); done
  done
  [ ${#_keep[@]} -gt 0 ] || die "CELLS_ONLY='$CELLS_ONLY' matches no cell (see: bash run_E3.sh plan)"
  CELLS=("${_keep[@]}")
fi
read -r -a SEED_ARR <<< "$SEEDS"
NCELL=${#CELLS[@]}
NSEED=${#SEED_ARR[@]}
NTASK=$((NCELL * NSEED))

if [ "${1:-}" != "check" ]; then
    for f in cpg_small.py cpg_e0.py cpg_e1.py cpg_e2.py cpg_e3.py; do
        [ -f "$HERE/$f" ] || die "$f must be in $HERE -- run: bash run_E3.sh check"
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
    for f in cpg_small.py cpg_e0.py cpg_e1.py cpg_e2.py cpg_e3.py run_E3.sh summarize_small.py; do
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
    if [ "$bad" = 0 ]; then echo; echo "  All good. Next: bash run_E3.sh smoke"
    else echo; echo "  Not ready: conda activate nest, and put the files listed above in $HERE"; return 1; fi
}

cmd_plan () {
    say "E3 grid: $NCELL cells x $NSEED seeds = $NTASK tasks per pass (x2 passes), ${SIM_S}s each"
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
    if [ -z "${E3_AWAKE:-}" ] && command -v caffeinate >/dev/null 2>&1; then
        echo "  (macOS: re-running under caffeinate -i so the Mac does not sleep mid-run)"
        E3_AWAKE=1 exec caffeinate -i bash "$HERE/run_E3.sh" --local "$which"
    fi
    echo "  Serial run, 30 s each: ~5-10 s wall per ladder run, ~20 s per control run."
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
    command -v sbatch >/dev/null || die "sbatch not found -- use: bash run_E3.sh --local"
    mkdir -p "$ROOT"
    local exp="ALL,SEEDS=$SEEDS,SIM_S=$SIM_S,ROOT=$ROOT,LAMBDA=$LAMBDA,PY=$PY"
    exp="$exp,STRICT_NEST=$STRICT_NEST,CELLS_ONLY=${CELLS_ONLY:-}"
    local ja jb
    ja=$(sbatch --parsable -p "$PARTITION" --time="$TIME" --array=0-$((NTASK - 1)) \
         --export="$exp" "$HERE/run_E3.sh" task A) || die "pass-A submit failed"
    jb=$(sbatch --parsable -p "$PARTITION" --time="$TIME" --array=0-$((NTASK - 1)) \
         --export="$exp" --dependency=afterany:"$ja" "$HERE/run_E3.sh" task B) \
         || die "pass-B submit failed"
    echo "  pass A: job $ja   pass B: job $jb   ->  when B is done: bash run_E3.sh summarize"
}

cmd_smoke () {
    local spec="n100_kmin3c_stdpoff|n100|off|indegree|0|--indeg-min 3 --indeg-min-conserve" x
    say "Smoke test: n100_kmin3c_stdpoff, 5 s, seed ${SEED_ARR[0]}, pass A then B"
    ROOT="$ROOT/smoke" SIM_S=5 FORCE=1
    run_one A "$spec" "${SEED_ARR[0]}" "$X0F_STOCK" "$X0F_STOCK" || die "pass A failed"
    x=$(prescribe n100_kmin3c_stdpoff) || die "prescription failed"
    run_one B "$spec" "${SEED_ARR[0]}" "${x% *}" "${x#* }" || die "pass B failed"
    "$PY" "$SUMMARY" --root "$ROOT" --ref-prefix ctrl_n480 --prefix "$ROOT/E3_summary" >/dev/null \
        && echo "  summarizer ok" || echo "  (summarizer FAILED)"
    echo "  ok -- outputs in $ROOT   (5 s, 1 seed: NOT evidence)"
}

cmd_summarize () {
    [ -f "$SUMMARY" ] || die "summarize_small.py not found: $SUMMARY"
    "$PY" "$SUMMARY" --root "$ROOT" --ref-prefix ctrl_n480 --prefix "$ROOT/E3_summary" \
        --t0 15 ${BASELINE_H5:+--baseline "$BASELINE_H5"}
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
    *) sed -n '13,21p' "${BASH_SOURCE[0]}"; exit 1 ;;
esac
