from pdbfixer import PDBFixer
from openmm.app import PDBFile, PDBxFile, Modeller, ForceField
from openff.toolkit import Molecule, Topology
from rdkit import Chem


import importlib 
import openff_agent.utils
importlib.reload(openff_agent.utils)


def find_target_product(
    products,
    resid_number: int = 126,
    tag_smiles: str = "C(=O)[O-]",
):
    """
    在 reaction 结果 products 里：
    - 对每个 rdmol 施加一个“肽键切割 + 端基封盖”反应
    - 在切出来的 capped 残基里找到：
        - 指定残基号 resid_number
        - 且 SMILES 中包含 tag_smiles 的那个 capped AA

    返回:
      (full_product_mol, capped_aa_mol, residue_number, chain_id)
    """
    rxn = Chem.rdChemReactions.ReactionFromSmarts(
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


def get_offmol(
    capped_aa: Chem.Mol,
    residue_number: int,
    chain_id: str,
) -> Molecule:
    """
    把 RDKit 的 capped AA 转成 OFFMol，并：
    - sanitize
    - generate conformers
    - 识别 ACE / NME cap，并在对应原子上打 residue metadata
    - 添加默认层级（chain / residue 等）
    """

    # 1. Sanitize RDKit mol
    Chem.SanitizeMol(capped_aa)

    # 2. RDKit Mol -> OpenFF Molecule
    offmol = Molecule.from_rdkit(
        capped_aa,
        allow_undefined_stereo=True,  # 你那边立体化学有坑，先放宽
    )

    # 3. 生成构象（OpenFF 需要 3D 坐标时会用）
    offmol.generate_conformers()

    # 4. 给 ACE / NME 端基打 residue metadata
    _annotate_caps_with_metadata(offmol, residue_number, chain_id)

    # 5. 添加默认层级信息（chain / residue 层次）
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


def make_ptm_smiles(noptm_pdb, ptm_smiles, smarts_rxn, output_folder):
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

    
def build_ptm_ff():
    pass


if __name__ == "__main__":

    noptm_pdb = "/eagle/projects/FoundEpidem/xlian/Agent/OpenffAgent/tmp/noptm/pred.model_idx_0.cif"
    ptm_smiles = 'CC(=O)C(=O)[O-]'
    smarts_rxn = Chem.rdChemReactions.ReactionFromSmarts(
        '[C:3][C:4][C:5][C:6][N;+1:1]([H:7])([H:8])[H:9].C[C:2](=O)C(=O)[O-]' # 给赖氨酸侧链的碳也编上号
        '>>'
        '[C:3][C:4][C:5][C:6][N;+0:1]=[C:2](C([H:7])([H:8])[H:9])(C(=O)[O-])'
    )
    output_folder = "../tmp/"

    make_ptm_smiles(noptm_pdb, ptm_smiles, smarts_rxn, output_folder)

