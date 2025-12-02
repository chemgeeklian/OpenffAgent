unfixed_pdb_path = '/eagle/projects/FoundEpidem/xlian/Agent/OpenffAgent/tmp/ptm/pred.model_idx_0.cif'
fixed_pdb_path = '/eagle/projects/FoundEpidem/xlian/Agent/OpenffAgent/tmp/ptm/pred.model_idx_0_fixed.pdb'

kpi_smiles_path = '/eagle/projects/FoundEpidem/xlian/Agent/OpenffAgent/tmp/noptm_H_add_KPI_smiles.txt'

from openff.toolkit import Molecule
from openmm.app import PDBFile
from rdkit import Chem
from pdbfixer import PDBFixer
import openff_agent.utils as utils

# get smiles from txt
smiles = open(kpi_smiles_path, 'r').read()

print('number of atoms in smiles: ')
print(len(Chem.AddHs(Chem.MolFromSmiles(smiles)).GetAtoms()))

utils.fix_pdb(unfixed_pdb_path, fixed_pdb_path, resid_to_rm_atom={'126': ['HO1', 'H2']})
pdb = PDBFile(fixed_pdb_path)

num_atoms = sum(1 for _ in pdb.topology.atoms())
print("number of atoms in PDB:", num_atoms)

mol_withKPI_from_smiles = Molecule.from_pdb_and_smiles(fixed_pdb_path, smiles, allow_undefined_stereo=True)
print(f"✅ Loaded molecule with from_pdb_and_smiles()")

# diff protein_with_KPI.pdb /eagle/projects/FoundEpidem/xlian/bio_ai_agent_genslm_example/energy_minim/openff_test/input/2v82_protein_with_H.pdb > diff.diff

# can load after removing H2 and HO1 from KPI