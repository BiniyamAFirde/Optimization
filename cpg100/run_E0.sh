#!/bin/bash -l
#SBATCH --job-name=cpg_E0
#SBATCH --output=slurm_E0_%A_%a.out
#SBATCH --error=slurm_E0_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=01:30:00
#SBATCH --mem=4G
#
# run_E0.sh -- roadmap node E0: baseline + control at 60 s, STDP off and on
# =========================================================================
#
#   bash run_E0.sh check             files + python packages + NEST present?  (run first)
#   bash run_E0.sh plan              cells, task count (no NEST)
#   bash run_E0.sh counts            build each network once, print ACTUAL counts
#   bash run_E0.sh smoke             5 s, ctrl_n480_stdpoff + base_k100_stdpon, A then B
#   bash run_E0.sh submit            Slurm: pass-A array, pass-B array (afterany)
#   bash run_E0.sh --local [A|B|all] the same thing serially (Mac)
#   bash run_E0.sh summarize         table + stability + force traces (pass B only)
#
#   (Slurm itself calls:  run_E0.sh task A|B   with SLURM_ARRAY_TASK_ID)
#
# ── GRID: 5 cells x 8 seeds x 2 passes = 80 runs, 60 s each ─────────────
#   cell                    preset     network                             STDP       chunk
#   base_k100_stdpoff       base_k100  k=1.00, 1200 neurons, inh-comp 1.0  off        100 ms
#   base_k100_stdpon        base_k100  "                                   on (1e-5)  100 ms
#   ctrl_n480_stdpoff       ctrl_n480  k=0.20, 480 neurons, inh-comp 3.75  off        100 ms
#   ctrl_n480_stdpon        ctrl_n480  "                                   on (1e-5)  100 ms
#   chk20_ctrl_n480_stdpoff ctrl_n480  "                                   off         20 ms
#
#   chunk = trace sample interval AND muscle/Ia update step. 100 ms reproduces the
#   model of record (plot1 samples every ~86.7 ms); cpg_small's default 20 ms does
#   not (1-seed probe: corr_RG -0.64 vs -0.95). chk20_* measures that at 8 seeds.
#   It is NOT part of the E0 verdict (its name does not match the targets).
#
# ── TWO PASSES ───────────────────────────────────────────────────────────
#   A  stock gate X0_F = 1.42 -> only the RG-F band is used
#   B  X0_F = median over seeds of gate_diag's X0F_rec (per cell)
#   K_F 30, K_E 15, X0_E 1.10 never change. Judge pass B only.
#
# Every run: --threads 1. Keep all 64 runs on ONE machine / ONE NEST build.
# Budget (1 thread, measured on a 2026 x86 node, NEST 3.10): ctrl_n480 60 s ~ 1-2 min,
# base_k100 60 s ~ 20 min (1200 neurons, 118k synapses, ~23 M spikes/s). --time is
# sized for base_k100 with 4x headroom (override: TIME=hh:mm:ss). cpg_e0.py empties
# the spike recorders after each read, so memory stays < 1 GB (it was ~20 GB).

set -uo pipefail

HERE="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
PY="${PY:-python3}"
DRIVER="${DRIVER:-$HERE/cpg_e0.py}"
GATE_DIAG="${GATE_DIAG:-$HERE/gate_diag.py}"
SUMMARY="${SUMMARY:-$HERE/summarize_small.py}"
ROOT="${ROOT:-$HERE/results_E0}"
SEEDS="${SEEDS:-12345 22345 32345 42345 52345 62345 72345 82345}"
SIM_S="${SIM_S:-60}"
LAMBDA="${LAMBDA:-1e-5}"
PARTITION="${PARTITION:-acc}"
TIME="${TIME:-01:30:00}"
FORCE="${FORCE:-0}"
STRICT_NEST="${STRICT_NEST:-0}"      # 1 = abort unless NEST is exactly 3.9.0
X0F_STOCK=1.42
# E0 reproduction targets (pass B mean corr_F, leg L), tolerance 0.02
TARGETS="${TARGETS:-base_k100=-0.953,ctrl_n480=-0.94}"

# name | preset | stdp | chunk_ms
CELLS=(
  "base_k100_stdpoff|base_k100|off|100"
  "base_k100_stdpon|base_k100|on|100"
  "ctrl_n480_stdpoff|ctrl_n480|off|100"
  "ctrl_n480_stdpon|ctrl_n480|on|100"
  "chk20_ctrl_n480_stdpoff|ctrl_n480|off|20"
)
say () { printf '\n\033[1m== %s\033[0m\n' "$*"; }
die () { printf '\n[FATAL] %s\n' "$*" >&2; exit 1; }

if [ -n "${CELLS_ONLY:-}" ]; then            # e.g. CELLS_ONLY="ctrl_n480_stdpoff"
  _keep=()
  for c in "${CELLS[@]}"; do
    for w in $CELLS_ONLY; do [ "${c%%|*}" = "$w" ] && _keep+=("$c"); done
  done
  [ ${#_keep[@]} -gt 0 ] || die "CELLS_ONLY='$CELLS_ONLY' matches no cell (see: bash run_E0.sh plan)"
  CELLS=("${_keep[@]}")
fi
read -r -a SEED_ARR <<< "$SEEDS"
NCELL=${#CELLS[@]}
NSEED=${#SEED_ARR[@]}
NTASK=$((NCELL * NSEED))

[ -f "$DRIVER" ] || die "driver not found: $DRIVER"
if [ "${1:-}" != "check" ]; then
    [ -f "$HERE/cpg_small.py" ] || die "cpg_small.py must sit next to cpg_e0.py ($HERE) -- run: bash run_E0.sh check"
fi
cell_field () { echo "$1" | cut -d'|' -f"$2"; }

# prescribe <cell> -> X0_F for pass B (median of pass-A X0F_rec over seeds)
prescribe () {
    local cell=$1
    "$PY" - "$ROOT/passA" "$cell" "$GATE_DIAG" "$SIM_S" <<'PY'
import glob, json, os, statistics, sys
d, cell, gd, sim_s = sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4])
transient = min(10.0, 0.2 * sim_s)
vals, src = [], ""
files = sorted(glob.glob(os.path.join(d, f"{cell}_seed*.h5")))
if os.path.isfile(gd):
    sys.path.insert(0, os.path.dirname(os.path.abspath(gd)))
    import gate_diag as G
    for p in files:
        try:
            meta, tr, _ = G.load(p, "L")
            vals.append(G.analyse(meta, tr, transient)["X0F_rec"])
        except Exception as e:
            print(f"[prescribe] skip {os.path.basename(p)}: {e!r}", file=sys.stderr)
    src = "gate_diag"
if not vals:
    for p in sorted(glob.glob(os.path.join(d, f"{cell}_seed*.json"))):
        if p.endswith(".counts.json"):
            continue
        try:
            vals.append(json.load(open(p))["rgf_band"]["X0F_rec"])
        except Exception:
            pass
    src = "pass-A json (same formula as gate_diag)"
if not vals:
    sys.exit(f"[prescribe] no pass-A results for {cell} in {d}")
x0 = statistics.median(vals)
print(f"[prescribe] {cell}: X0_F = {x0:.3f}  (median of {len(vals)} seeds, {src})", file=sys.stderr)
print(f"{x0:.3f}")
PY
}

# run_one <A|B> <cell-spec> <seed> <x0f>
run_one () {
    local pass=$1 spec=$2 seed=$3 x0f=$4
    local cell preset stdp chunk
    cell=$(cell_field "$spec" 1); preset=$(cell_field "$spec" 2); stdp=$(cell_field "$spec" 3)
    chunk=$(cell_field "$spec" 4)
    local outdir="$ROOT/pass$pass" tag="${cell}_seed${seed}"
    mkdir -p "$outdir"
    if [ -s "$outdir/$tag.json" ] && [ "$FORCE" != "1" ]; then
        echo "  [skip] pass $pass $tag"; return 0
    fi
    local extra=()
    [ "$STRICT_NEST" = "1" ] && extra+=(--strict-nest)
    echo "  [run ] pass $pass $tag  preset=$preset stdp=$stdp chunk=${chunk}ms X0_F=$x0f  (${SIM_S}s)"
    "$PY" -u "$DRIVER" --preset "$preset" --out "$outdir" --tag "$tag" --cell "$cell" \
        --pass-label "$pass" --stdp "$stdp" --stdp-lambda "$LAMBDA" \
        --sim-s "$SIM_S" --seed "$seed" --x0f "$x0f" --x0e 1.10 --threads 1 \
        --chunk-ms "$chunk" ${extra[@]+"${extra[@]}"} > "$outdir/$tag.log" 2>&1
    local rc=$?
    if [ $rc -ne 0 ]; then
        echo "        FAILED rc=$rc -- tail of $outdir/$tag.log:"
        tail -5 "$outdir/$tag.log" | sed 's/^/        | /'
        return $rc
    fi
    grep -E '^\[(done|E0)\]|^\[counts\] (Izh|syn)' "$outdir/$tag.log" | sed 's/^/        /'
}

task () {
    local pass=$1 idx=$2
    [ "$idx" -lt "$NTASK" ] || die "task index $idx >= $NTASK"
    local spec="${CELLS[$((idx / NSEED))]}" seed="${SEED_ARR[$((idx % NSEED))]}"
    local cell; cell=$(cell_field "$spec" 1)
    local x0f=$X0F_STOCK
    if [ "$pass" = "B" ]; then
        x0f=$(prescribe "$cell") || die "no pass-A prescription for $cell"
        mkdir -p "$ROOT/passB"; echo "$x0f" > "$ROOT/passB/$cell.x0f"
    fi
    run_one "$pass" "$spec" "$seed" "$x0f"
}

cmd_check () {
    say "Pre-flight: $(uname -s) $(uname -m), bash $BASH_VERSION, python = $(command -v "$PY")"
    local bad=0 f
    for f in cpg_small.py cpg_e0.py run_E0.sh summarize_small.py gate_diag.py; do
        if [ -f "$HERE/$f" ]; then
            printf '  ok       %s\n' "$f"
        elif [ "$f" = gate_diag.py ]; then
            printf '  missing  %s  (optional: pass B then uses the same formula from the pass-A JSON)\n' "$f"
        else
            printf '  MISSING  %s\n' "$f"; bad=1
        fi
    done
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
    v = nest.__version__
    print(f"  ok       nest {v}" + ("" if v == "3.9.0" else
          "   <-- not 3.9.0: fine for a self-consistent local set, but do not mix with FINDINGS/cluster numbers"))
    nest.ResetKernel(); n = nest.Create("izhikevich", 5); r = nest.Create("spike_recorder")
    nest.Connect(n, r); n.set(I_e=10.0); nest.Simulate(100.0)
    print(f"  ok       NEST smoke simulation ({r.get('n_events')} spikes)")
except Exception as e:
    ok = False; print(f"  MISSING  nest  ({e.__class__.__name__}: {e})")
sys.exit(0 if ok else 1)
PY
    if [ "$bad" = 0 ]; then
        echo; echo "  All good. Next: bash run_E0.sh smoke"
    else
        echo
        echo "  Not ready. Usual fixes:"
        echo "    - activate the NEST env first:   conda activate nest     (the prompt must NOT say (base))"
        echo "      (no env yet?  bash mac_nest.sh setup   -- or:"
        echo "       conda create -y -n nest -c conda-forge nest-simulator=3.9.0 python=3.11 numpy h5py matplotlib)"
        echo "    - copy cpg_small.py and gate_diag.py from the project into: $HERE"
        return 1
    fi
}

cmd_plan () {
    say "E0 grid: $NCELL cells x $NSEED seeds = $NTASK tasks per pass (x2 passes), ${SIM_S}s each"
    local c; for c in "${CELLS[@]}"; do
        printf '  %-24s preset=%-10s stdp=%-4s chunk=%s ms\n' "$(cell_field "$c" 1)" \
            "$(cell_field "$c" 2)" "$(cell_field "$c" 3)" "$(cell_field "$c" 4)"
    done
    echo "  lambda=$LAMBDA  seeds: $SEEDS"
    echo "  targets: $TARGETS (pass B, leg L, tol 0.02)"
    echo "  results -> $ROOT/passA, $ROOT/passB"
}

cmd_counts () {
    say "Actual counts NEST builds (build only, seed ${SEED_ARR[0]})"
    mkdir -p "$ROOT/counts"
    local p; for p in base_k100 ctrl_n480; do
        "$PY" "$DRIVER" --preset "$p" --build-only --out "$ROOT/counts" --tag "$p" \
            --seed "${SEED_ARR[0]}" --threads 1 > "$ROOT/counts/$p.build.log" 2>&1 \
            || { echo "  $p: BUILD FAILED -- last lines of $ROOT/counts/$p.build.log:";
                 tail -4 "$ROOT/counts/$p.build.log" | sed 's/^/     | /'; continue; }
        grep -E '^\[counts\] (Izh|syn)|^!!' "$ROOT/counts/$p.build.log" | sed "s/^/  $p: /"
    done
}

cmd_local () {
    local which=${1:-all} i
    if [ -z "${E0_AWAKE:-}" ] && command -v caffeinate >/dev/null 2>&1; then
        echo "  (macOS: re-running under caffeinate -i so the Mac does not sleep mid-run)"
        E0_AWAKE=1 exec caffeinate -i bash "$HERE/run_E0.sh" --local "$which"
    fi
    echo "  Serial run. Budget per 60 s run, 1 thread: ctrl ~1-2 min, base_k100 ~15-25 min."
    echo "  Finished runs are skipped on restart, so Ctrl-C is safe; re-run the same command to resume."
    if [ "$which" = "A" ] || [ "$which" = "all" ]; then
        say "Pass A (stock gate X0_F=$X0F_STOCK) -- $NTASK runs, serial"
        for ((i = 0; i < NTASK; i++)); do task A "$i"; done
    fi
    if [ "$which" = "B" ] || [ "$which" = "all" ]; then
        say "Pass B (each cell at its own prescribed X0_F) -- $NTASK runs, serial"
        for ((i = 0; i < NTASK; i++)); do task B "$i"; done
    fi
    [ "$which" = "all" ] && cmd_summarize
}

cmd_submit () {
    command -v sbatch >/dev/null || die "sbatch not found -- use: bash run_E0.sh --local"
    mkdir -p "$ROOT"
    local exp="ALL,SEEDS=$SEEDS,SIM_S=$SIM_S,ROOT=$ROOT,LAMBDA=$LAMBDA,PY=$PY"
    exp="$exp,STRICT_NEST=$STRICT_NEST,CELLS_ONLY=${CELLS_ONLY:-}"
    local ja jb
    ja=$(sbatch --parsable -p "$PARTITION" --time="$TIME" --array=0-$((NTASK - 1)) \
         --export="$exp" "$HERE/run_E0.sh" task A) || die "pass-A submit failed"
    jb=$(sbatch --parsable -p "$PARTITION" --time="$TIME" --array=0-$((NTASK - 1)) \
         --export="$exp" --dependency=afterany:"$ja" "$HERE/run_E0.sh" task B) \
         || die "pass-B submit failed"
    echo "  pass A: job $ja  ($NTASK tasks)"
    echo "  pass B: job $jb  (starts when every pass-A task has ended)"
    echo "  watch:  squeue -u \$USER     when B is done:  bash run_E0.sh summarize"
}

cmd_smoke () {
    say "Smoke test: ctrl_n480_stdpoff + base_k100_stdpon, 5 s, seed ${SEED_ARR[0]}, pass A then B"
    ROOT="$ROOT/smoke" SIM_S=5 FORCE=1
    local spec cell x0
    for spec in "ctrl_n480_stdpoff|ctrl_n480|off|100" "base_k100_stdpon|base_k100|on|100"; do
        cell=$(cell_field "$spec" 1)
        run_one A "$spec" "${SEED_ARR[0]}" "$X0F_STOCK" || die "pass A failed ($cell)"
        x0=$(prescribe "$cell") || die "prescription failed ($cell)"
        run_one B "$spec" "${SEED_ARR[0]}" "$x0" || die "pass B failed ($cell)"
    done
    "$PY" "$SUMMARY" --root "$ROOT" --targets "$TARGETS" --ref-prefix base_k100 \
        --prefix "$ROOT/E0_summary" >/dev/null \
        && echo "  summarizer ok -> $ROOT/E0_summary.csv" || echo "  (summarizer FAILED)"
    echo "  ok -- outputs in $ROOT   (5 s, 1 seed: NOT evidence)"
}

cmd_summarize () {
    [ -f "$SUMMARY" ] || die "summarize_small.py not found: $SUMMARY"
    "$PY" "$SUMMARY" --root "$ROOT" --targets "$TARGETS" --ref-prefix base_k100 \
        --prefix "$ROOT/E0_summary" ${BASELINE_H5:+--baseline "$BASELINE_H5"}
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
    *) sed -n '13,24p' "${BASH_SOURCE[0]}"; exit 1 ;;
esac
