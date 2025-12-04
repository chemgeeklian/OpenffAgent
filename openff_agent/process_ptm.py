'''
Generate force field for modified amino acids.
This only needs to run once to get FF for the modified AA...

Steps:
1. From predicted PTM structure, apply reaction to get products.
2. Find target product.
3. Generate force field for target product.
'''

import openff_agent.utils
import openff

from rdkit import Chem
from rdkit.Chem import rdChemReactions
from chemper.mol_toolkits import mol_toolkit
from chemper.graphs.single_graph import SingleGraph
from openff.toolkit import Molecule, Topology, ForceField
from openff.units import unit

def find_target_product(
    products,
    resid_number: int = 126,
    tag_smiles: str = "C(=O)[O-]",
):
    """
    In the reaction products:
    - apply a "peptide bond cleavage + terminal capping" reaction to each rdmol
    - find the capped AA that:
        - has the specified residue number resid_number
        - and has the tag_smiles in its SMILES

    Returns:
      (full_product_mol, capped_aa_mol, residue_number, chain_id)
    """
    rxn = rdChemReactions.ReactionFromSmarts(
        "[C](=[O])[N:3]([H:4])[C:5][C:6](=[O:7])[N]"
        " >> "
        "[CH3]C(=O)[N:3]([H:4])[C:5][C:6](=[O:7])N([H])[CH3]"
    )

    for index, rdmol in enumerate(products):
        cut_products = rxn.RunReactants([rdmol[0]])
        for product in cut_products:
            capped_aa = product[0]
            for atom in capped_aa.GetAtoms():
                resinfo = atom.GetPDBResidueInfo()
                if resinfo is None:
                    continue

                resnum = resinfo.GetResidueNumber()
                chain_id = resinfo.GetChainId()
                resname = resinfo.GetResidueName()

                if resnum != resid_number:
                    continue

                smiles = Chem.MolToSmiles(capped_aa)
                if tag_smiles in smiles:
                    # 找到了目标残基
                    print(f"Hit at index={index}, res={resname}{resnum}, smiles={smiles}")
                    return rdmol[0], capped_aa, resnum, chain_id

    raise RuntimeError(
        f"No capped AA found with resid_number={resid_number} and tag_smiles='{tag_smiles}'"
    )


def get_capped_ncaa_offmol(
    capped_aa: Chem.Mol,
    residue_number: int,
    chain_id: str,
) -> Molecule:
    """
    translate above
    Convert the RDKit capped AA to an OpenFF Molecule, and:
    - sanitize
    - generate conformers
    - identify ACE / NME caps and tag them with residue metadata
    - add default hierarchy (chain / residue etc.)
    """

    # 1. Sanitize RDKit mol
    Chem.SanitizeMol(capped_aa)

    # 2. RDKit Mol -> OpenFF Molecule
    offmol = Molecule.from_rdkit(
        capped_aa,
        allow_undefined_stereo=True,
    )

    offmol.generate_conformers()

    # Assign residue metadata for ACE / NME caps
    _annotate_caps_with_metadata(offmol, residue_number, chain_id)

    offmol.add_default_hierarchy_schemes()

    return offmol


def _annotate_caps_with_metadata(
    offmol: Molecule,
    residue_number: int,
    chain_id: str,
):
    """
    在 offmol 里用 SMARTS 找 ACE / NME，
    给对应匹配到的原子补充 residue metadata：
      - residue_name: ACE / NME
      - residue_number: residue_number ± 1
      - insertion_code: " "
      - chain_id: chain_id
    """
    cap_definitions = [
        # N 端 ACE: CH3-CO-N
        ("ACE", "[C:1]([H:2])([H:3])([H:4])[C:5](=[O:6])N", -1),
        # C 端 NME: -CONH-CH3
        ("NME", "[C:1]([H:2])([H:3])([H:4])[N:5]([H:6])C", +1),
    ]

    for cap_resname, smarts, resoffset in cap_definitions:
        matches = offmol.chemical_environment_matches(smarts, unique=True)

        for match in matches:
            for idx in match:
                atom = offmol.atom(idx)
                # 只在还没有 residue_name 的情况下覆盖，避免污染真正的主链残基
                if atom.metadata.get("residue_name") is None:
                    atom.metadata["residue_name"] = cap_resname
                    atom.metadata["residue_number"] = residue_number + resoffset
                    atom.metadata["chain_id"] = chain_id
                    atom.metadata["insertion_code"] = " "


def make_ptm_products(noptm_pdb, 
                      ptm_smiles: str, 
                      smarts_rxn: rdChemReactions.ChemicalReaction):
    if '.cif' in noptm_pdb:
        noptm_pred_pdb_fixed = noptm_pdb.replace(".cif", "_fixed.pdb")
    elif '.pdb' in noptm_pdb:
        noptm_pred_pdb_fixed = noptm_pdb.replace(".pdb", "_fixed.pdb")
    else:
        raise ValueError("Input structure must be .pdb or .cif format")

    openff_agent.utils.fix_pdb(noptm_pdb, noptm_pred_pdb_fixed)

    topology = Topology.from_pdb(noptm_pred_pdb_fixed)
    molecule = topology.molecule(0)

    modification_mal = Molecule.from_smiles(ptm_smiles)

    protein = molecule.to_rdkit()
    modification = modification_mal.to_rdkit()   # reagent

    products = smarts_rxn.RunReactants([protein, modification])

    print("Number of products:", len(products))

    return products

    
def build_ptm_ff(ncaa_capped: Molecule,
                 output_path):
    ncaa_capped.assign_partial_charges('am1bcc') # conda install -y -c conda-forge ambertools !!!
    # sage_ff14sb = ForceField('openff-2.0.0.offxml', 'ff14sb_off_impropers_0.0.3.offxml')
    sage_ff14sb = ForceField("openff-2.2.0.offxml", "ff14sb_off_impropers_0.0.4.offxml")

    ncaa_uncapped_indices = []
    for res in ncaa_capped.residues:
        if res.residue_name not in ["ACE", "NME"]:
            ncaa_uncapped_indices += [atom.molecule_atom_index for atom in res.atoms]

    formal_charge = 0. * unit.elementary_charge
    total_partial_charges = 0. * unit.elementary_charge

    for i in ncaa_uncapped_indices:
        atom = ncaa_capped.atoms[i]
        formal_charge += atom.formal_charge
        total_partial_charges += atom.partial_charge
        
    charge_offset = (formal_charge - total_partial_charges) / len(ncaa_uncapped_indices)

    mol = mol_toolkit.Mol(ncaa_capped.to_rdkit())
    graph = SingleGraph(mol, ncaa_uncapped_indices, layers=1)
    charges_smirks = graph.as_smirks()

    charges = [None for i in ncaa_uncapped_indices]

    # Chemper's `SingleGraph` class produces smirks index `:1` corresponding to `ncaa_uncapped_indices[0]`,
    # so we can iterate over the indices, put the appropriate charge in place, and apply the offset we 
    # calculated earlier
    for smirks_atom_index, molecule_atom_index in enumerate(ncaa_uncapped_indices, start=1):
        charges[smirks_atom_index-1] = ncaa_capped.partial_charges[molecule_atom_index] + charge_offset
        
    # Make sure we haven't doubly assigned any charges
    assert None not in charges

    # Delete any charges we might've put in the last time we ran this cell
    try:
        del sage_ff14sb['LibraryCharges'].parameters[charges_smirks]
    except openff.toolkit.utils.exceptions.ParameterLookupError:
        pass

    # Add the new charges to the library with the ChemPer SMIRKS
    sage_ff14sb['LibraryCharges'].add_parameter({'smirks': charges_smirks,  'charge': charges})
    sage_ff14sb.to_file(output_path)


if __name__ == "__main__":

    noptm_pdb = "../tmp/noptm/pred.model_idx_0.cif"
    ptm_smiles = 'CC(=O)C(=O)[O-]'
    smarts_rxn = rdChemReactions.ReactionFromSmarts(
        '[C:3][C:4][C:5][C:6][N;+1:1]([H:7])([H:8])[H:9].C[C:2](=O)C(=O)[O-]'
        '>>'
        '[C:3][C:4][C:5][C:6][N;+0:1]=[C:2](C([H:7])([H:8])[H:9])(C(=O)[O-])'
    )
    output_folder = "../output/"

    products = make_ptm_products(noptm_pdb, ptm_smiles, smarts_rxn)

    full_product_rdmol, capped_aa_rdmol, resid_number, chain_id = find_target_product(
        products,
        resid_number=126,
        tag_smiles="C(=O)[O-]",
    )

    ncaa_capped = get_capped_ncaa_offmol(capped_aa_rdmol, resid_number, chain_id)
    build_ptm_ff(ncaa_capped, output_folder + "KPI.offxml")

    print("Finished generating PTM force field.")