#!/usr/bin/env python3
import argparse
import gemmi
from pdbfixer import PDBFixer

def cif2pdb(cif_file: str, pdb_file: str=None):

    if pdb_file is None:
        pdb_file = cif_file.replace(".cif", ".pdb")

    structure = gemmi.read_structure(cif_file)
    structure.write_pdb(pdb_file)
    print(f"Converted {cif_file} → {pdb_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("cif_file", type=str, help="input CIF file")
    parser.add_argument("pdb_file", type=str, nargs="?", help="output PDB file (optional)")
    args = parser.parse_args()

    cif_file = args.cif_file
    pdb_file = args.pdb_file if args.pdb_file else cif_file.replace(".cif", ".pdb")

    cif2pdb(cif_file, pdb_file)