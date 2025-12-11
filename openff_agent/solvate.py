from openff.units import Quantity, unit
from openmm.app import PDBFile
from rdkit import Chem

import numpy
import openmm
import openmm.app
import openmm.unit
from openff.toolkit import ForceField, Molecule, Topology
from openff.units import Quantity, unit
from openff.units.openmm import ensure_quantity
from openmm.app import Topology as OpenMMTopology

from openff.interchange import Interchange

OPENMM_IONS = {
    "Li+": "[#3+1]",
    "Na+": "[#11+1]",
    "K+": "[#19+1]",
    "Rb+": "[#37+1]",
    "Cs+": "[#55+1]",
    "F-": "[#9-1]",
    "Cl-": "[#17-1]",
    "Br-": "[#35-1]",
    "I-": "[#53-1]",
}

def solvate_topology(
    topology: Topology,
    method: str = "pdbfixer",
    box_vectors: Quantity | None = Quantity(5.0 * numpy.ones(3), unit.nanometer),
    **kwargs,
) -> Topology:
    if method in ["pdbfixer", "openmm"]:
        boxSize = openmm.unit.Quantity(openmm.Vec3(*box_vectors.m_as(unit.nanometer)), openmm.unit.nanometer)

        if method == "pdbfixer":
            openmm_topology, openmm_positions = _solvate_pdbfixer(
                topology.to_openmm(),
                topology.get_positions().to_openmm(),
                boxSize=boxSize,
                **kwargs,
            )
        else:
            openmm_topology, openmm_positions = _solvate_openmm(
                topology.to_openmm(),
                topology.get_positions().to_openmm(),
                boxSize=boxSize,
                **kwargs,
            )

        unique_molecules: list[Molecule] = [*topology.unique_molecules]
        unique_molecules.append(Molecule.from_mapped_smiles("[H:2][O:1][H:3]"))

        if "positiveIon" in kwargs:
            unique_molecules.append(Molecule.from_smiles(OPENMM_IONS[kwargs["positiveIon"]]))

        if "negativeIon" in kwargs:
            unique_molecules.append(Molecule.from_smiles(OPENMM_IONS[kwargs["negativeIon"]]))

        new_topology = Topology.from_openmm(
            openmm_topology,
            unique_molecules=unique_molecules,
        )

        new_topology.set_positions(ensure_quantity(openmm_positions, "openff"))

        return new_topology


def _solvate_pdbfixer(
    topology: OpenMMTopology,
    positions: openmm.unit.Quantity,
    **kwargs,
) -> tuple[OpenMMTopology, openmm.unit.Quantity]:
    """
    Add solvent and ions using PDBFixer.

    https://htmlpreview.github.io/?https://github.com/openmm/pdbfixer/blob/master/Manual.html

    """
    import pdbfixer
    import os

    with open("_tmp.pdb", "w") as _file:
        openmm.app.PDBFile.writeFile(topology, positions, _file)

    pdb_object = pdbfixer.PDBFixer("_tmp.pdb")
    pdb_object.addSolvent(**kwargs)
    os.remove("_tmp.pdb")

    modeller = openmm.app.Modeller(pdb_object.topology, pdb_object.positions)

    # remove 3 CL-, gross!
    cl_residues = [r for r in modeller.topology.atoms() if r.name == 'Cl']
    modeller.delete(cl_residues[:3])

    return modeller.topology, modeller.positions


def _solvate_openmm(
    topology: OpenMMTopology,
    positions: openmm.unit.Quantity,
    box_vectors: openmm.unit.Quantity,
    forcefield: openmm.app.ForceField | None = None,
    **kwargs,
) -> tuple[OpenMMTopology, openmm.unit.Quantity]:
    if not forcefield:
        import pdbfixer

        forcefield = pdbfixer.PDBFixer._createForceField(topology)

    modeller = openmm.app.Modeller(topology, positions)
    modeller.addSolvent(
        forcefield,
        **kwargs,
    )


def insert_molecule_and_remove_clashes(
    topology: Topology,
    insert: Molecule,
    radius: Quantity = 1.5 * unit.angstrom,
    keep: list[Molecule] = [],
) -> Topology:
    """
    Add a molecule to a copy of the topology, removing any clashing molecules.

    The molecule will be added to the end of the topology. A new topology is
    returned; the input topology will not be altered. All molecules that
    clash will be removed, and each removed molecule will be printed to stdout.
    Users are responsible for ensuring that no important molecules have been
    removed; the clash radius may be modified accordingly.

    Parameters
    ==========
    top
        The topology to insert a molecule into
    insert
        The molecule to insert
    radius
        Any atom within this distance of any atom in the insert is considered
        clashing.
    keep
        Keep copies of these molecules, even if they're clashing
    """
    # We'll collect the molecules for the output topology into a list
    new_top_mols = []
    # A molecule's positions in a topology are stored as its zeroth conformer
    insert_coordinates = insert.conformers[0][:, None, :]
    for molecule in topology.molecules:
        if any(keep_mol.is_isomorphic_with(molecule) for keep_mol in keep):
            new_top_mols.append(molecule)
            continue

        ''' 
        molecule_coordinates = molecule.conformers[0][None, :, :]
       
        diff_matrix = molecule_coordinates - insert_coordinates

        # np.linalg.norm doesn't work on Pint quantities 😢
        working_unit = unit.nanometer
        distance_matrix = np.linalg.norm(diff_matrix.m_as(working_unit), axis=-1) * working_unit

        if distance_matrix.min() > radius:
            # This molecule is not clashing, so add it to the topology
            new_top_mols.append(molecule)
        else:
            print(f"Removed {molecule.to_smiles()} molecule")
        '''

        new_top_mols.append(molecule)

    # Insert the ligand at the end
    new_top_mols.append(ligand)

    # This pattern of assembling a topology from a list of molecules
    # ends up being much more efficient than adding each molecule
    # to a new topology one at a time
    new_top = Topology.from_molecules(new_top_mols)

    # Don't forget the box vectors!
    new_top.box_vectors = topology.box_vectors
    return new_top


TMP_DIR = "../tmp"   # <-- 单独抽出来

ligand_path = f"{TMP_DIR}/ptm_lig/e4p_docked.sdf"
pdb_path = f"{TMP_DIR}/ptm_lig/pred.model_idx_0_fixed_chainA.pdb"
kpi_smiles_path = f"{TMP_DIR}/ptm_product.smiles"

# Load a molecule from a SDF file
ligand = Molecule.from_file(ligand_path)

# Print out a SMILES code for the ligand
print(ligand.to_smiles(explicit_hydrogens=False))

# get smiles from txt
smiles = open(kpi_smiles_path, 'r').read()

print('number of atoms in smiles: ')
print(len(Chem.AddHs(Chem.MolFromSmiles(smiles)).GetAtoms()))

pdb = PDBFile(pdb_path)
num_atoms = sum(1 for _ in pdb.topology.atoms())
print("number of atoms in PDB:", num_atoms)

mol_withKPI_from_smiles = Molecule.from_pdb_and_smiles(pdb_path, smiles, allow_undefined_stereo=True)
print(f"✅ Loaded molecule with from_pdb_and_smiles()")

topology = mol_withKPI_from_smiles.to_topology()
print('number of molecules in topology: ', topology.n_molecules)


solvated_topology = solvate_topology(
    topology,
    box_vectors=Quantity(7.0 * numpy.ones(3), unit.nanometer),
    positiveIon="Na+",
    negativeIon="Cl-",
    ionicStrength=0.15 * openmm.unit.molar,
)


top = insert_molecule_and_remove_clashes(solvated_topology, ligand)
n_na = sum(1 for atom in top.atoms if atom.symbol.upper() == "NA")
n_cl = sum(1 for atom in top.atoms if atom.symbol.upper() == "CL")

print("Na+", n_na)
print("Cl-", n_cl)

with open(f"{TMP_DIR}/solvated_topology.json", "w") as f:
    print(top.to_json(), file=f)

print('wrote to json.')

'''
sage_ff14sb = ForceField( "../output/KPI.offxml")
interchange = sage_ff14sb.create_interchange(top)

import time
# Length of the simulation.
num_steps = 1000  # number of integration steps to run

# Logging options.
trj_freq = 10  # number of steps per written trajectory frame
data_freq = 10  # number of steps per written simulation statistics

# Integration options
time_step = 2 * openmm.unit.femtoseconds  # simulation timestep
temperature = 300 * openmm.unit.kelvin  # simulation temperature
friction = 1 / openmm.unit.picosecond  # friction constant

integrator = openmm.LangevinMiddleIntegrator(temperature, friction, time_step)

simulation = interchange.to_openmm_simulation(integrator=integrator)

openmm_system = interchange.to_openmm()
openmm_topology = interchange.to_openmm_topology()
openmm_positions = interchange.positions.to_openmm()

simulation.context.setVelocitiesToTemperature(temperature)
simulation.context.setPositions(openmm_positions)
simulation.minimizeEnergy()

state = simulation.context.getState(getPositions=True)
min_positions = state.getPositions()

# Write PDB
with open("minimized.pdb", "w") as f:
    PDBFile.writeFile(openmm_topology, min_positions, f)
'''