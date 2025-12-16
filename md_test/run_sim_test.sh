#!/bin/bash -l
# Run `sim_test.py` multiple times for two systems and post-process distances.

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

run_one() {
	local src_topology=$1
	local dest_run_dir=$2

	mkdir -p "$dest_run_dir"
	# sim_test expects a topology file at TMP_DIR/solvated_topology.json
	cp "$src_topology" "$dest_run_dir/solvated_topology.json"

	echo "Running sim_test.py with tmp_dir=$dest_run_dir"
	python sim_test.py --tmp_dir "$dest_run_dir" &> "$dest_run_dir/sim_test.log"
}

# Run repeats sequentially (safe and predictable)
for i in $(seq 1 $REPEATS); do
	echo "== WT repetition $i =="
	run_one "$WT_SRC" "$WT_BASE/run_$i"
done

for i in $(seq 1 $REPEATS); do
	echo "== EC03_1 repetition $i =="
	run_one "$EC03_SRC" "$EC03_BASE/run_$i"
done

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
