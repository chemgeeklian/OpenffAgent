#!/usr/bin/env python

import os, sys
from pathlib import Path

import parsl
from parsl import python_app

# --------------------------------------------------
# Inputs
# --------------------------------------------------
run_dir = Path('/nfs/ml_lab/projects/ml_lab/xlian/IDEAL/nibrin/nbn_chai_16samples/runs/')
pdbs = [
    Path('/nfs/ml_lab/projects/ml_lab/xlian/IDEAL/nibrin/nbn_chai_16samples/selected_pdb/gpu1_pred_model_idx_2.pdb'),
    Path('/nfs/ml_lab/projects/ml_lab/xlian/IDEAL/nibrin/nbn_chai_16samples/selected_pdb/gpu2_pred_model_idx_3.pdb')
]
paths = [run_dir / pdb.stem for pdb in pdbs]

# Safety checks
if not run_dir.is_dir():
    raise NotADirectoryError(f"Run directory not found: {run_dir}")

for pdb in pdbs:
    if not pdb.is_file():
        raise FileNotFoundError(f"PDB file not found: {pdb}")

from parsl.config import Config
from parsl.executors.threads import ThreadPoolExecutor

config = Config(
            executors=[
                        ThreadPoolExecutor(
                                        label='local_threads',
                                                    max_threads=4   # 允许最多 4 个并行任务
                                                            )
                            ]
            )

parsl.load(config)   # local parallelism
# If you want to use more threads:
# from parsl.configs.local_threads import config
# parsl.load(config)

# --------------------------------------------------
# Simulation App
# --------------------------------------------------
@python_app
def run_simulation(path, pdb, device_ids):
    import sys
    sys.path.insert(0,
        '/eagle/projects/FoundEpidem/xlian/IDEAL/nibrin/openmm_mat/molecular-dynamics/src/'
    )
    from molecular_simulations.simulate.omm_simulator import Simulator

    print(f"Starting simulation for {path} on GPUs {device_ids}")

    sim = Simulator(
        path,
        device_ids=device_ids,
        temperature=310,
        equil_steps=1_250_000,
        prod_steps=500_000_000,
        eq_reporter_frequency=25_000,
        prod_reporter_frequency=50_000
    )
    sim.run()
    return f"Simulation completed for {path}"

# --------------------------------------------------
# Launch two simulations simultaneously
# --------------------------------------------------
device_assignments = [[0,1], [2,3]]

futures = [
    run_simulation(path, pdb, device_ids)
    for path, pdb, device_ids in zip(paths, pdbs, device_assignments)
]

# Wait for both
results = [f.result() for f in futures]
print(results)
print("All simulations completed.")

