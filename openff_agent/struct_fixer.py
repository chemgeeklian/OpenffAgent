ptm_pdb_path = '../tmp/ptm_lig/pred.model_idx_0.cif'
noptm_pdb_path = "../tmp/noptm/pred.model_idx_0.cif"

from openff.toolkit import Molecule
from openmm.app import PDBFile
from rdkit import Chem
from rdkit.Chem import rdChemReactions
import openff_agent.utils as utils
from openff_agent.process_ptm import make_ptm_products, find_target_product


fixed_pdb_path = ptm_pdb_path.replace(".cif", "_fixed.pdb")
splited_enzyme_pdb_path = fixed_pdb_path.replace(".pdb", f"_chainA.pdb")

utils.fix_pdb(ptm_pdb_path, fixed_pdb_path, resid_to_rm_atom={'126': ['HO1', 'H2']})
utils.split_pdb(fixed_pdb_path)

pdb = PDBFile(splited_enzyme_pdb_path)

num_atoms = sum(1 for _ in pdb.topology.atoms())
print("number of atoms in PDB:", num_atoms)


ptm_smiles = 'CC(=O)C(=O)[O-]'
smarts_rxn = rdChemReactions.ReactionFromSmarts(
    '[C:3][C:4][C:5][C:6][N;+1:1]([H:7])([H:8])[H:9].C[C:2](=O)C(=O)[O-]'
    '>>'
    '[C:3][C:4][C:5][C:6][N;+0:1]=[C:2](C([H:7])([H:8])[H:9])(C(=O)[O-])'
)

products = make_ptm_products(noptm_pdb_path, ptm_smiles, smarts_rxn)

full_product_rdmol, capped_aa_rdmol, resid_number, chain_id = find_target_product(
    products,
    resid_number=126,
    tag_smiles="C(=O)[O-]",
)

smiles = Chem.MolToSmiles(full_product_rdmol)

print('number of atoms in smiles: ')
print(len(Chem.AddHs(Chem.MolFromSmiles(smiles)).GetAtoms()))


mol_withKPI_from_smiles = Molecule.from_pdb_and_smiles(splited_enzyme_pdb_path, 
                                                       smiles, 
                                                       allow_undefined_stereo=True)
print(f"✅ Loaded molecule with from_pdb_and_smiles()")


# can load after removing H2 and HO1 from KPI