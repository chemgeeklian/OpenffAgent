pdb_path = '/eagle/projects/FoundEpidem/xlian/Agent/OpenffAgent/tmp/ptm_H.pdb'
kpi_smiles_path = '/eagle/projects/FoundEpidem/xlian/Agent/OpenffAgent/tmp/predmodel_KPI_smiles.txt'

from openff.toolkit import Molecule, Topology
from openmm.app import PDBFile, Modeller, ForceField as OpenMMForceField
from pdbfixer import PDBFixer
import rdkit
from rdkit import Chem
import os

# get smiles from txt
smiles = open(kpi_smiles_path, 'r').read()

print('number of atoms in smiles: ')
print(len(Chem.AddHs(Chem.MolFromSmiles(smiles)).GetAtoms()))

pdb = PDBFile(pdb_path)
num_atoms = sum(1 for _ in pdb.topology.atoms())
print("number of atoms in PDB:", num_atoms)

mol_withKPI_from_smiles = Molecule.from_pdb_and_smiles(pdb_path, smiles, allow_undefined_stereo=True)
print(f"✅ Loaded molecule with from_pdb_and_smiles()")

# diff protein_with_KPI.pdb /eagle/projects/FoundEpidem/xlian/bio_ai_agent_genslm_example/energy_minim/openff_test/input/2v82_protein_with_H.pdb > diff.diff