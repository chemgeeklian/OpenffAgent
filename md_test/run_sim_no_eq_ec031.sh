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
        python sim_no_eq.py --tmp_dir "$out" &> "$out/sim_no_eq.log"
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