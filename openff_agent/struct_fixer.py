unfixed_pdb_path = '/eagle/projects/FoundEpidem/xlian/Agent/OpenffAgent/tmp/ptm_lig/pred.model_idx_0.cif'
kpi_smiles_path = '/eagle/projects/FoundEpidem/xlian/Agent/OpenffAgent/tmp/noptm_H_add_KPI_smiles.txt'

from openff.toolkit import Molecule
from openmm.app import PDBFile
from rdkit import Chem
import openff_agent.utils as utils


fixed_pdb_path = unfixed_pdb_path.replace(".cif", "_fixed.pdb")
splited_enzyme_pdb_path = fixed_pdb_path.replace(".pdb", f"_chainA.pdb")

utils.fix_pdb(unfixed_pdb_path, fixed_pdb_path, resid_to_rm_atom={'126': ['HO1', 'H2']})
utils.split_pdb(fixed_pdb_path)

pdb = PDBFile(splited_enzyme_pdb_path)

num_atoms = sum(1 for _ in pdb.topology.atoms())
print("number of atoms in PDB:", num_atoms)


# get smiles from txt
# TODO: implement codes to do rxn and get SMILES

smiles = open(kpi_smiles_path, 'r').read()

print('number of atoms in smiles: ')
print(len(Chem.AddHs(Chem.MolFromSmiles(smiles)).GetAtoms()))

mol_withKPI_from_smiles = Molecule.from_pdb_and_smiles(splited_enzyme_pdb_path, 
                                                       smiles, 
                                                       allow_undefined_stereo=True)
print(f"✅ Loaded molecule with from_pdb_and_smiles()")


# can load after removing H2 and HO1 from KPI