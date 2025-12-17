from rdkit import Chem
from rdkit.Chem import AllChem

def dock(target_lig_coord_pdb: str, smiles: str, output_sdf: str):
    """
    Superpose ligand generated from SMILES onto an existing ligand pose stored in a PDB file.
    target_lig_coord_pdb: PDB containing ONLY the ligand coordinates in binding pose.
    """
    print(f"Loading reference pose: {target_lig_coord_pdb}")
    ref = Chem.MolFromPDBFile(target_lig_coord_pdb, sanitize=False, removeHs=False)
    if ref is None:
        raise ValueError("Failed to read reference ligand PDB.")

    ref = Chem.AddHs(ref, addCoords=True)

    print(f"Generating new ligand from SMILES: {smiles}")
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")

    mol = Chem.AddHs(mol)

    # Generate 3D conformer
    AllChem.EmbedMolecule(mol, randomSeed=123)
    AllChem.MMFFOptimizeMolecule(mol)

    # ------------------------------------------
    # Atom mapping (first try exact substructure)
    # ------------------------------------------
    match = mol.GetSubstructMatch(ref)

    # fallback: map by element symbol in order
    if not match:
        print("[WARN] No substructure match; fallback to element-wise matching.")
        match = []
        used = set()
        for atom_ref in ref.GetAtoms():
            sym = atom_ref.GetSymbol()
            for idx, atom_new in enumerate(mol.GetAtoms()):
                if idx in used:
                    continue
                if atom_new.GetSymbol() == sym:
                    match.append(idx)
                    used.add(idx)
                    break

    match = tuple(match)
    ref_idx = tuple(range(len(match)))

    # ------------------------------------------
    # Alignment
    # ------------------------------------------
    print("Aligning new ligand to reference pose...")
    AllChem.AlignMol(mol, ref, atomMap=list(zip(match, ref_idx)))

    # ------------------------------------------
    # Write output
    # ------------------------------------------
    w = Chem.SDWriter(output_sdf)
    w.write(mol)
    w.close()
    print(f"Docked ligand written to: {output_sdf}")


if __name__ == "__main__":
    TMP_DIR = "../wt_tmp"
    
    target_lig_coord_pdb = f"{TMP_DIR}/ptm_lig/pred.model_idx_0_fixed_chainB.pdb"
    smiles = "O=C[C@H](O)[C@H](O)COP(=O)([O-])[O-]"
    output_sdf = f"{TMP_DIR}/ptm_lig/e4p_docked.sdf"
    dock(target_lig_coord_pdb, smiles, output_sdf)