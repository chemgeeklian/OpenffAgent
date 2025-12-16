import mdtraj as md
import numpy as np

'''
This script recenters a molecular dynamics trajectory on the protein's center of mass (COM).
It performs the following steps:
1) Loads the trajectory and topology.
2) Images the molecules to account for periodic boundary conditions.
3) Recenters each frame on the protein's COM.
4) Saves the modified trajectory and a PDB of the first frame.
'''

pdb = "../tmp/minimized.pdb"
dcd = "../md_test/run/eq.dcd"
output_dcd_path = "../md_test/run/eq_centered.dcd"

print("Loading trajectory...")
traj = md.load(dcd, top=pdb)

# Select protein atoms
protein_idx = traj.topology.select("protein")
if protein_idx.size == 0:
    raise RuntimeError("No protein atoms found (topology issue?).")

print(f"Protein atoms: {protein_idx.size}")

# 1) initial imaging
print("Imaging molecules...")
traj.image_molecules(inplace=True)

# 2) recenter to protein COM
print("Centering on protein COM...")
for i, frame in enumerate(traj):
    xyz = frame.xyz[0]  # shape = (n_atoms, 3)
    
    protein_xyz = xyz[protein_idx]
    com = protein_xyz.mean(axis=0)

    box = frame.unitcell_lengths[0]  # (Lx, Ly, Lz)
    box_center = box / 2.0

    shift = box_center - com
    frame.xyz[0] = xyz + shift

# 3) wrap again
print("Final wrapping...")
traj.image_molecules(inplace=True)

# save
traj.save("../md_test/run/eq_centered.dcd")
traj[0].save_pdb("../md_test/run/eq_frame0_centered.pdb")

print("Saved eq_centered.dcd and eq_frame0_centered.pdb")
