#!/bin/bash -l
#PBS -l select=1:system=polaris
#PBS -l walltime=1:00:00
#PBS -q debug
#PBS -A FoundEpidem
#PBS -l filesystems=home:eagle
#PBS -N offsim

set -euo pipefail

module load conda
conda activate openff312
cd /home/xlian/OpenffAgent/md_test || exit 1

# ============================================================
# CONFIG
# ============================================================
REPEATS=4
ROOT=/home/xlian/OpenffAgent

# ============================================================
# 🔥 SYSTEMS TO RUN (ONLY TOUCH THIS BLOCK)
# ============================================================
SYSTEMS=(
  "WT:$ROOT/wt_tmp/solvated_topology.json:$ROOT/wt_runs"
  "EC03_1:$ROOT/tmp/solvated_topology.json:$ROOT/ec03_1_runs"
)

SYSTEMS=(
  "EC03_1:$ROOT/tmp/solvated_topology.json:$ROOT/ec03_1_runs"
)

# ============================================================
# GPU DETECTION
# ============================================================
detect_gpus() {
    if [ -n "${CUDA_VISIBLE_DEVICES:-}" ]; then
        echo "$CUDA_VISIBLE_DEVICES" | tr ',' ' '
    elif command -v nvidia-smi >/dev/null 2>&1; then
        nvidia-smi --query-gpu=index --format=csv,noheader,nounits
    else
        echo 0
    fi
}

GPUS=($(detect_gpus))
NGPU=${#GPUS[@]}
echo "Using GPUs: ${GPUS[*]} (NGPU=$NGPU)"

# ============================================================
# RUNNER
# ============================================================
run_one() {
    local src=$1
    local out=$2
    local gpu=$3

    mkdir -p "$out"
    cp "$src" "$out/solvated_topology.json"

    echo "→ $out on GPU $gpu"
    CUDA_VISIBLE_DEVICES=$gpu \
        python sim_test.py --tmp_dir "$out" &> "$out/sim_test.log"
}

# ============================================================
# RUN SIMULATIONS
# ============================================================
job=0
for entry in "${SYSTEMS[@]}"; do
    IFS=: read -r NAME SRC BASE <<< "$entry"
    echo "Starting $NAME runs"

    for i in $(seq 1 $REPEATS); do
        gpu=${GPUS[$((job % NGPU))]}
        run_one "$SRC" "$BASE/run_$i" "$gpu" &
        job=$((job+1))
        (( job % NGPU == 0 )) && wait
    done
    wait
done

# ============================================================
# POST-PROCESS
# ============================================================
PLOTS_DIR="$ROOT/plots"
mkdir -p "$PLOTS_DIR"

python - <<'PY'
from pathlib import Path
from compute_distance import process_systems, plot_results

ROOT = Path('/home/xlian/OpenffAgent')
REPEATS = 4

systems = {}

# MUST match bash SYSTEMS
system_defs = [
    # ("WT", "wt_runs"),
    ("EC03_1", "ec03_1_runs"),
]

for name, base in system_defs:
    for i in range(1, REPEATS + 1):
        systems[f"{name}_run{i}"] = ROOT / base / f"run_{i}" / "run"

print("Processing:")
for k, v in systems.items():
    print(" ", k, "->", v)

results = process_systems(systems)
out = ROOT / "plots/LIG_C1x_to_KPI_C1_distance_vs_time_ns.png"
plot_results(results, frame_dt_ps=10.0, filename=str(out), show=False)
print("Saved:", out)
PY