#!/usr/bin/env python3
"""
Production-only MD (NO equilibration).
- Minimize
- Initialize velocities
- Run production
"""

import argparse
import time
from pathlib import Path

from openmm.app import *
from openmm import *
from openmm.unit import *

from openff.toolkit import ForceField, Topology


# ============================================================================
# ARGS
# ============================================================================
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--tmp_dir", type=str, required=True,
        help="Path to temporary working directory."
    )
    return parser.parse_args()


args = parse_args()
TMP_DIR = Path(args.tmp_dir)

# ============================================================================
# PATHS
# ============================================================================
topology_json = TMP_DIR / "solvated_topology.json"
forcefield_file = "../output/KPI.offxml"

output_dir = TMP_DIR / "run"
output_dir.mkdir(exist_ok=True)

prod_log = output_dir / "prod.log"
prod_dcd = output_dir / "prod.dcd"
prod_chkpt = output_dir / "prod.chk"
prod_state = output_dir / "prod.state"
min_pdb = output_dir / "minimized.pdb"

# ============================================================================
# SIMULATION PARAMETERS
# ============================================================================
temperature = 300.0 * kelvin
friction = 1.0 / picosecond

prod_dt = 0.001 * picoseconds
prod_steps = 10_000_000
prod_freq = 10_000

# ============================================================================
# LOAD TOPOLOGY
# ============================================================================
print("Loading OpenFF topology...")

top2load = Topology.from_json(topology_json.read_text())
ff = ForceField(forcefield_file)
interchange = ff.create_interchange(top2load)

openmm_topology = interchange.to_openmm_topology()
openmm_positions = interchange.positions.to_openmm()
openmm_system = interchange.to_openmm()

print(f"System has {openmm_topology.getNumAtoms()} atoms")

# ============================================================================
# SETUP SYSTEM (PRODUCTION SYSTEM DIRECTLY)
# ============================================================================
print("Setting up production system (no restraints, NPT)...")

openmm_system.addForce(
    MonteCarloBarostat(1.0 * bar, temperature)
)

integrator = LangevinMiddleIntegrator(
    temperature,
    friction,
    prod_dt
)

platform = Platform.getPlatformByName("CUDA")

simulation = Simulation(
    openmm_topology,
    openmm_system,
    integrator,
    platform
)

# ============================================================================
# MINIMIZATION
# ============================================================================
print("Minimizing energy...")
simulation.context.setPositions(openmm_positions)

simulation.minimizeEnergy()

state = simulation.context.getState(getPositions=True, getEnergy=True)
print("Minimized potential energy:", state.getPotentialEnergy())

with open(min_pdb, "w") as f:
    PDBFile.writeFile(openmm_topology, state.getPositions(), f)

# ============================================================================
# INITIALIZE VELOCITIES
# ============================================================================
print(f"Initializing velocities at {temperature}...")
simulation.context.setVelocitiesToTemperature(temperature)

# ============================================================================
# REPORTERS
# ============================================================================
simulation.reporters.extend([
    DCDReporter(str(prod_dcd), prod_freq),
    StateDataReporter(
        str(prod_log),
        prod_freq,
        step=True,
        potentialEnergy=True,
        temperature=True,
        volume=True,
        speed=True,
        progress=True,
        remainingTime=True,
        totalSteps=prod_steps,
        separator="\t",
    ),
    CheckpointReporter(str(prod_chkpt), prod_freq * 10),
])

# ============================================================================
# PRODUCTION
# ============================================================================
print(
    f"Running production MD: "
    f"{prod_steps} steps "
    f"({prod_steps * prod_dt.value_in_unit(picoseconds) / 1000:.1f} ns)"
)

start = time.time()
simulation.step(prod_steps)
end = time.time()

simulation.saveState(str(prod_state))

elapsed = end - start
print("Production complete")
print(f"Elapsed: {elapsed/60:.1f} min")
print(f"Performance: {prod_steps * prod_dt.value_in_unit(picoseconds) / elapsed:.2f} ns/day")

print("\nOUTPUT:")
print("  minimized:", min_pdb)
print("  prod:", prod_dcd, prod_log, prod_chkpt, prod_state)
