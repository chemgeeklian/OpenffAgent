#!/usr/bin/env python

import os, sys
from pathlib import Path
import parsl
from parsl import python_app

# --------------------------------------------------
# Inputs
# --------------------------------------------------
run_dir = Path('run/')
paths = [run_dir]

# Safety checks
if not run_dir.is_dir():
    raise NotADirectoryError(f"Run directory not found: {run_dir}")

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
def run_simulation(path, device_ids=None):
    import sys
    from molecular_simulations.simulate.omm_simulator import Simulator

    sim = Simulator(
        path,
        platform='OpenCL',
        temperature=300,
        equil_steps=1_250_000,
        prod_steps=25_000_000, # 100ns
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
    run_simulation(path) for path in paths
]

# Wait for both
results = [f.result() for f in futures]
print(results)
print("All simulations completed.")

