import os
from pathlib import Path
from openff.toolkit import Molecule, Topology, ForceField
from openmm.app import PDBFile, PDBxFile, Modeller, ForceField
import os
from pdbfixer import PDBFixer


# ---------------------------
# 1. Load structure (CIF or PDB)
# ---------------------------
def load_structure(file_path: Path | str) -> PDBFixer:
    file_path = Path(file_path)
    suffix = file_path.suffix.lower()

    if suffix == ".pdb":
        fixer = PDBFixer(filename=str(file_path))

    elif suffix in (".cif", ".mmcif"):
        cif = PDBxFile(str(file_path))
        fixer = PDBFixer(topology=cif.topology, positions=cif.positions)

    else:
        raise ValueError("Unsupported file format: must be .pdb or .cif")

    return fixer


# ---------------------------
# 2. Fix missing residues & atoms
# ---------------------------
def fix_structure(fixer: PDBFixer) -> PDBFixer:
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    return fixer


# ---------------------------
# 3. Add hydrogens based on pH
# ---------------------------
def add_hydrogens(fixer: PDBFixer, ph: float = 7.4):
    modeller = Modeller(fixer.topology, fixer.positions)
    forcefield = ForceField("amber14-all.xml", "amber14/tip3pfb.xml")

    modeller.addHydrogens(forcefield, pH=ph)
    return modeller.topology, modeller.positions


# ---------------------------
# 4. Save structure as PDB or CIF
# ---------------------------
def save_structure(topology, positions, out_path: Path | str):
    out_path = Path(out_path)

    if out_path.suffix.lower() == ".pdb":
        with open(out_path, "w") as f:
            PDBFile.writeFile(topology, positions, f, keepIds=True)

    elif out_path.suffix.lower() == ".cif":
        with open(out_path, "w") as f:
            PDBxFile.writeFile(topology, positions, f)

    else:
        raise ValueError("Output must be .pdb or .cif")

    print(f"💾 Saved: {out_path}")


# ---------------------------
# 5. Main orchestrator (decoupled)
# ---------------------------
def prepare_structure(
    in_file: Path | str,
    out_file: Path | str,
    ph: float = 7.4,
):
    """
    Load → fix → protonate → save
    """

    fixer = load_structure(in_file)
    fixer = fix_structure(fixer)

    topology, positions = add_hydrogens(fixer, ph=ph)
    save_structure(topology, positions, out_file)

    return Path(out_file)
