import argparse
from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout
import time
from copy import deepcopy
from pathlib import Path
from openff.toolkit import ForceField, Topology


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tmp_dir", type=str, required=True,
                        help="Path to temporary working directory.")
    return parser.parse_args()

args = parse_args()
TMP_DIR = Path(args.tmp_dir)

topology_json = TMP_DIR / "solvated_topology.json"
forcefield_file = "../output/KPI.offxml"
output_dir = TMP_DIR / "run"
output_dir.mkdir(exist_ok=True)

# Simulation parameters
temperature = 300  # Kelvin
friction = 1.0  # 1/ps

# Equilibration parameters
equil_steps = 1_000_000  # Total equilibration steps (1 ns at 1 fs)
equil_dt = 0.001  # ps (1 fs)
equil_freq = 100_000  # Reporting frequency
force_constant = 10.0  # kcal/mol/Å^2 for restraints
equil_cycles = 3  # Number of NPT cycles after restraint relaxation

# Production parameters
prod_steps = 10_000_000  # Production steps (10 ns at 1 fs)
prod_dt = 0.001  # ps (1 fs)
prod_freq = 10_000  # Reporting frequency

# Output files
eq_log = output_dir / "eq.log"
eq_dcd = output_dir / "eq.dcd"
eq_chkpt = output_dir / "eq.chk"
eq_state = output_dir / "eq.state"
prod_log = output_dir / "prod.log"
prod_dcd = output_dir / "prod.dcd"
prod_chkpt = output_dir / "prod.chk"
prod_state = output_dir / "prod.state"

# ============================================================================
# LOAD TOPOLOGY AND CREATE INTERCHANGE
# ============================================================================
print("Loading OpenFF topology and creating Interchange...")
with open(topology_json, "r") as f:
    json_str = f.read()

top2load = Topology.from_json(json_str)
sage_ff14sb = ForceField(forcefield_file)
interchange = sage_ff14sb.create_interchange(top2load)

# Get OpenMM objects from Interchange
openmm_topology = interchange.to_openmm_topology()
openmm_positions = interchange.positions.to_openmm()
openmm_system = interchange.to_openmm()

print(f"System has {openmm_topology.getNumAtoms()} atoms")

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================
def get_backbone_indices(topology):
    """Get indices of protein backbone atoms (CA, C, N, O)."""
    indices = []
    for atom in topology.atoms():
        if atom.residue.name in ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU',
                                 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE',
                                 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'KPI', 'UNK']:
            if atom.name in ['CA', 'C', 'N', 'O']:
                indices.append(atom.index)
    print(f"Found {len(indices)} backbone atoms for restraints")
    return indices

def add_backbone_restraints(system, positions, topology, indices, k=10.0):
    """Add harmonic position restraints to backbone atoms."""
    force = CustomExternalForce("k*periodicdistance(x, y, z, x0, y0, z0)^2")

    force_amount = k * kilocalories_per_mole / angstroms**2
    force.addGlobalParameter("k", force_amount)
    force.addPerParticleParameter("x0")
    force.addPerParticleParameter("y0")
    force.addPerParticleParameter("z0")

    atoms_list = list(topology.atoms())
    for i, (atom_crd, atom) in enumerate(zip(positions, atoms_list)):
        if atom.index in indices:
            force.addParticle(i, atom_crd.value_in_unit(nanometers))

    restrained_system = deepcopy(system)
    restrained_system.addForce(force)

    return restrained_system

def heating_protocol(simulation, integrator, target_temp=300, start_temp=5):
    """Gradually heat system from start_temp to target_temp."""
    print(f"\nHeating from {start_temp}K to {target_temp}K...")
    simulation.context.setVelocitiesToTemperature(start_temp * kelvin)

    integrator.setTemperature(start_temp * kelvin)
    mdsteps = 100_000
    heat_steps = 1_000
    length = mdsteps // heat_steps
    tstep = (target_temp - start_temp) / length

    for i in range(length):
        simulation.step(heat_steps)
        temp = start_temp + tstep * (1 + i)
        if temp > target_temp:
            temp = target_temp
        integrator.setTemperature(temp * kelvin)

        if (i + 1) % 10 == 0:
            print(f"  Step {(i+1)*heat_steps}/{mdsteps}, T = {temp:.1f}K")

    print(f"Heating complete at {target_temp}K")
    return simulation, integrator

def restraint_relaxation(simulation, k_initial=10.0, n_levels=5, steps_per_level=None):
    """Gradually relax restraints in NVT ensemble."""
    print(f"\nRelaxing restraints from {k_initial} kcal/mol/Å^2 to 0...")
    simulation.context.reinitialize(True)

    d_k = k_initial / n_levels

    for i in range(n_levels):
        simulation.step(steps_per_level)
        k = float(k_initial - (i * d_k))
        simulation.context.setParameter('k', k * kilocalories_per_mole / angstroms**2)
        print(f"  Level {i+1}/{n_levels}: k = {k:.2f} kcal/mol/Å^2, ran {steps_per_level} steps")

    # Final unrestrained NVT
    simulation.context.setParameter('k', 0)
    simulation.step(steps_per_level)
    print(f"  Unrestrained NVT: k = 0, ran {steps_per_level} steps")

    return simulation


# ============================================================================
# EQUILIBRATION
# ============================================================================
print("\n" + "="*70)
print("EQUILIBRATION PHASE")
print("="*70)

# Get backbone indices for restraints
backbone_indices = get_backbone_indices(openmm_topology)

# Add restraints to system
print("\nAdding backbone restraints...")
restrained_system = add_backbone_restraints(
    openmm_system,
    openmm_positions,
    openmm_topology,
    backbone_indices,
    k=force_constant
)

# Create integrator and simulation for equilibration
print("\nSetting up equilibration simulation...")
eq_integrator = LangevinMiddleIntegrator(
    temperature * kelvin,
    friction / picosecond,
    equil_dt * picoseconds
)

platform = Platform.getPlatformByName('CUDA')

eq_simulation = Simulation(
    openmm_topology,
    restrained_system,
    eq_integrator,
    platform=platform
)

# Set positions and minimize energy
print("\nMinimizing energy...")
eq_simulation.context.setPositions(openmm_positions)
initial_state = eq_simulation.context.getState(getEnergy=True, getPositions=True)
print(f"Initial potential energy: {initial_state.getPotentialEnergy()}")

eq_simulation.minimizeEnergy()

minimized_state = eq_simulation.context.getState(getEnergy=True, getPositions=True)
print(f"Minimized potential energy: {minimized_state.getPotentialEnergy()}")

# save minimized
with open(f"{output_dir}/minimized.pdb", "w") as f:
    PDBFile.writeFile(openmm_topology, minimized_state.getPositions(), f)

# Add reporters
print("\nAdding reporters...")
eq_simulation.reporters.append(
    StateDataReporter(
        str(eq_log),
        equil_freq,
        step=True,
        potentialEnergy=True,
        temperature=True,
        speed=True,
        progress=True,
        remainingTime=True,
        totalSteps=equil_steps + 100_000,  # heating + equilibration
        separator='\t'
    )
)
eq_simulation.reporters.append(DCDReporter(str(eq_dcd), equil_freq))

# Heating protocol
eq_simulation, eq_integrator = heating_protocol(
    eq_simulation,
    eq_integrator,
    target_temp=temperature,
    start_temp=5
)

# Restraint relaxation in NVT
steps_per_level = equil_steps // 5
eq_simulation = restraint_relaxation(
    eq_simulation,
    k_initial=force_constant,
    n_levels=5,
    steps_per_level=steps_per_level
)

# Add barostat for NPT equilibration
print("\nAdding barostat for NPT equilibration...")
barostat = MonteCarloBarostat(1.0 * bar, temperature * kelvin)
eq_simulation.system.addForce(barostat)

# Run NPT equilibration cycles
print(f"\nRunning {equil_cycles} NPT equilibration cycles...")
for cycle in range(equil_cycles):
    print(f"  NPT cycle {cycle+1}/{equil_cycles}")
    eq_simulation.step(steps_per_level)

# Save equilibration state and checkpoint
print("\nSaving equilibration state and checkpoint...")
eq_simulation.saveState(str(eq_state))
eq_simulation.saveCheckpoint(str(eq_chkpt))

print("\nEquilibration complete!")

# ============================================================================
# PRODUCTION
# ============================================================================
print("\n" + "="*70)
print("PRODUCTION PHASE")
print("="*70)

# Create new system with barostat (no restraints)
print("\nSetting up production simulation...")
prod_system = interchange.to_openmm()
prod_system.addForce(MonteCarloBarostat(1.0 * bar, temperature * kelvin))

# Create integrator and simulation for production
prod_integrator = LangevinMiddleIntegrator(
    temperature * kelvin,
    friction / picosecond,
    prod_dt * picoseconds
)

prod_simulation = Simulation(
    openmm_topology,
    prod_system,
    prod_integrator,
    platform=platform
)

# Load checkpoint from equilibration
print("\nLoading equilibration checkpoint...")
prod_simulation.loadCheckpoint(str(eq_chkpt))
state = prod_simulation.context.getState(getVelocities=True, getPositions=True)
positions = state.getPositions()
velocities = state.getVelocities()

prod_simulation.context.setPositions(positions)
prod_simulation.context.setVelocities(velocities)

# Add reporters
print("\nAdding reporters...")
prod_simulation.reporters.extend([
    DCDReporter(str(prod_dcd), prod_freq),
    StateDataReporter(
        str(prod_log),
        prod_freq,
        step=True,
        potentialEnergy=True,
        temperature=True,
        progress=True,
        remainingTime=True,
        speed=True,
        volume=True,
        totalSteps=prod_steps,
        separator='\t'
    ),
    CheckpointReporter(str(prod_chkpt), prod_freq * 10)
])

# Run production MD
print(f"\nRunning production MD for {prod_steps} steps ({prod_steps * prod_dt / 1000:.1f} ns)...")
start_time = time.time()
prod_simulation.step(prod_steps)
end_time = time.time()

# Save final state
print("\nSaving final state...")
prod_simulation.saveState(str(prod_state))

elapsed_time = end_time - start_time
print(f"\nProduction complete!")
print(f"Total time: {elapsed_time/60:.1f} minutes")
print(f"Performance: {prod_steps * prod_dt / elapsed_time:.2f} ns/day")

print("\n" + "="*70)
print("SIMULATION COMPLETE")
print("="*70)
print(f"\nOutput files:")
print(f"  Equilibration: {eq_log}, {eq_dcd}, {eq_chkpt}, {eq_state}")
print(f"  Production: {prod_log}, {prod_dcd}, {prod_chkpt}, {prod_state}")


