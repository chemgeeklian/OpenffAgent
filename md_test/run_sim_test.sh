#!/bin/bash -l
#PBS -l select=1:system=polaris
#PBS -l walltime=4:00:00
#PBS -q preemptable
#PBS -A FoundEpidem
#PBS -l filesystems=home:eagle
#PBS -N offsim

module load conda

conda activate openff312
cd /home/xlian/OpenffAgent/md_test || exit 1

# Configuration
REPEATS=3
ROOT=/home/xlian/OpenffAgent

# Source topology locations (existing prepared directories)
WT_SRC="$ROOT/wt_tmp/solvated_topology.json"
EC03_SRC="$ROOT/tmp/solvated_topology.json"

# Where to place repeated runs
WT_BASE="$ROOT/wt_runs"
EC03_BASE="$ROOT/ec03_1_runs"

mkdir -p "$WT_BASE" "$EC03_BASE"

echo "Starting simulation repeats: $REPEATS per system"

# Detect GPUs: use CUDA_VISIBLE_DEVICES if set, otherwise query nvidia-smi
detect_gpus() {
    local arr=()
    if [ -n "$CUDA_VISIBLE_DEVICES" ]; then
        IFS=',' read -ra tmp <<< "$CUDA_VISIBLE_DEVICES"
        for t in "${tmp[@]}"; do
            # trim whitespace
            t=$(echo "$t" | xargs)
            [ -n "$t" ] && arr+=("$t")
        done
    else
        if command -v nvidia-smi >/dev/null 2>&1; then
            while IFS= read -r idx; do
                arr+=("$idx")
            done < <(nvidia-smi --query-gpu=index --format=csv,noheader,nounits 2>/dev/null)
        fi
    fi
    if [ ${#arr[@]} -eq 0 ]; then
        echo "No GPUs detected; falling back to GPU 0 (may run on CPU)"
        arr=(0)
    fi
    echo "${arr[@]}"
}

GPUS=($(detect_gpus))
NGPU=${#GPUS[@]}
echo "Using GPUs: ${GPUS[*]} (NGPU=$NGPU)"

run_one() {
    local src_topology=$1
    local dest_run_dir=$2
    local gpu_id=$3

    mkdir -p "$dest_run_dir"
    cp "$src_topology" "$dest_run_dir/solvated_topology.json"

    echo "Running sim_test.py in $dest_run_dir on GPU $gpu_id"
    CUDA_VISIBLE_DEVICES=$gpu_id \
    python sim_test.py --tmp_dir "$dest_run_dir" &> "$dest_run_dir/sim_test.log"
}

echo "Starting WT runs"
job=0
for i in $(seq 1 $REPEATS); do
    gpu=${GPUS[$((job % NGPU))]}
    run_one "$WT_SRC" "$WT_BASE/run_$i" "$gpu" &
    job=$((job+1))
    (( job % NGPU == 0 )) && wait
done
wait

echo "Starting EC03_1 runs"
job=0
for i in $(seq 1 $REPEATS); do
    gpu=${GPUS[$((job % NGPU))]}
    run_one "$EC03_SRC" "$EC03_BASE/run_$i" "$gpu" &
    job=$((job+1))
    (( job % NGPU == 0 )) && wait
done
wait

echo "All simulations finished. Running distance post-processing..."

# Ensure plots directory exists
PLOTS_DIR="$ROOT/plots"
mkdir -p "$PLOTS_DIR"

# Post-process distances using compute_distance functions
python - <<'PY'
from pathlib import Path
from compute_distance import process_systems, plot_results

ROOT = Path('/home/xlian/OpenffAgent')
PLOTS = ROOT / 'plots'
systems = {}

REPEATS = 3
for i in range(1, REPEATS+1):
    systems[f'WT_run{i}'] = ROOT / f'wt_runs/run_{i}/run'
    systems[f'EC03_1_run{i}'] = ROOT / f'ec03_1_runs/run_{i}/run'

print('Processing distance for runs:')
for k, v in systems.items():
    print('  ', k, '->', v)

results = process_systems(systems)
out_file = PLOTS / 'LIG_C1x_to_KPI_C1_distance_vs_time_ns.png'
plot_results(results, frame_dt_ps=10.0, filename=str(out_file), show=False)
print('Distance processing and plotting complete. Saved to', out_file)
PY

echo "Done. Plots saved to $PLOTS_DIR"
