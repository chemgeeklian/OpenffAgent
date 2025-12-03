#!/usr/bin/env python3
import argparse
import gemmi
from pdbfixer import PDBFixer
from openmm.app import PDBFile, Modeller


def cif2pdb(cif_file: str, pdb_file: str=None):

    if pdb_file is None:
        pdb_file = cif_file.replace(".cif", ".pdb")

    structure = gemmi.read_structure(cif_file)
    structure.write_pdb(pdb_file)
    print(f"Converted {cif_file} → {pdb_file}")


def get_resid_id(fixer: PDBFixer, res_name: str="KPI"):
    '''Get the residue index of KPI in PDBFixer object.'''
    for residue in fixer.topology.residues():
        if residue.name == res_name:
            return residue.index
    return None


def _fix_resid(fixer: PDBFixer, resid_idx: int, rm_atoms: list):
    '''Remove specified atoms from a residue in PDBFixer object.'''
    for atom in fixer.topology.atoms():
        if atom.residue.index == resid_idx and atom.name in rm_atoms:
            fixer.topology.removeAtom(atom.index)


def fix_pdb(pdb_file: str, output_file: str=None, resid_to_rm_atom: dict=None, chain_id: list=["A"]):
    '''Use PDBFixer to fix common issues in PDB file.'''

    if output_file is None:
        output_file = pdb_file.replace(".pdb", "_fixed.pdb")

    fixer = PDBFixer(filename=pdb_file)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.4)

    modeller = Modeller(fixer.topology, fixer.positions)

    if resid_to_rm_atom:
        atoms_to_delete = []
        for chain in fixer.topology.chains():
            if chain.id in chain_id:
                for residue in chain.residues():
                    if str(residue.id) in resid_to_rm_atom:
                        for atom in residue.atoms():
                            if atom.name in resid_to_rm_atom[str(residue.id)]:
                                atoms_to_delete.append(atom)
        modeller.delete(atoms_to_delete)

    with open(output_file, 'w') as f:
        PDBFile.writeFile(modeller.topology, modeller.positions, f, keepIds=True)
    print(f"Fixed PDB file {pdb_file} → {output_file}")


def split_pdb(pdb_file:str, ):
    # split pdb by chains
    pdb = PDBFile(pdb_file)
    for chain in pdb.topology.chains():
        chain_id = chain.id
        output_file = pdb_file.replace(".pdb", f"_chain{chain_id}.pdb")
        atoms_to_keep = [
            atom for atom in pdb.topology.atoms()
            if atom.residue.chain.id == chain_id
        ]
        modeller = Modeller(pdb.topology, pdb.positions)
        modeller.delete([atom for atom in pdb.topology.atoms() if atom not in atoms_to_keep])
        with open(output_file, 'w') as f:
            PDBFile.writeFile(modeller.topology, modeller.positions, f, keepIds=True)
        print(f"Extracted chain {chain_id} to {output_file}")


if __name__ == "__main__":
    split_pdb('/eagle/projects/FoundEpidem/xlian/Agent/OpenffAgent/tmp/ptm_lig/pred.model_idx_0_fixed.pdb')