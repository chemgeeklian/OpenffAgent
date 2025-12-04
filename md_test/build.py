#!/usr/bin/env python

import os, sys
sys.path.insert(0, os.path.abspath('/eagle/projects/FoundEpidem/xlian/IDEAL/nibrin/openmm_mat/molecular-dynamics/src/'))
from molecular_simulations import PolarisSettings
import parsl
from parsl import python_app
#from glob import glob

yaml_file = 'polaris_debug.yaml'
run_dir = '/home/xlian/nbn_chai_16samples/runs/'
pdbs = ['/home/xlian/nbn_chai_16samples/selected_pdb/gpu1_pred_model_idx_2.pdb', '/home/xlian/nbn_chai_16samples/selected_pdb/gpu2_pred_model_idx_3.pdb']
paths = [run_dir + os.path.basename(pdb).split('.pdb')[-2] for pdb in pdbs]

# --------------------------------------------------
# Build Step
# --------------------------------------------------
#@python_app
def build_structure(path, pdb):
    import os, sys
    sys.path.insert(0, os.path.abspath('/eagle/projects/FoundEpidem/xlian/IDEAL/nibrin/openmm_mat/molecular-dynamics/src/'))
    from molecular_simulations.build.build_amber import ExplicitSolvent

    builder = ExplicitSolvent(path, pdb, protein=True, padding=20.)
    builder.build()
    print(f'Build completed for {path}')

# --------------------------------------------------
# Main Execution
# --------------------------------------------------

# Load Parsl config
# Load job settings from a YAML file, which specifies parameters 
# such as queue, walltime, number of workers, GPUs, etc.
# Create a Parsl Config object using these settings, with executors like HighThroughputExecutor, 
# WorkQueueExecutor, or ThreadPoolExecutor.
# Load the configuration into Parsl's runtime to initialize a pool of workers across nodes 
# (on Polaris) for executing the defined functions.

#settings = PolarisSettings.from_yaml(yaml_file)
#config = settings.config_factory(run_dir)
#parsl.load(config)

# Step 1: Build
build_futures = [build_structure(path, pdb) for path, pdb in zip(paths, pdbs)]
#_ = [f.result() for f in build_futures]
print("All structures built.")

# conda activate /eagle/projects/FoundEpidem/msinclair/conda_envs/mdsim/