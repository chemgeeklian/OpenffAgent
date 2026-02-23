import mdtraj as md
import numpy as np
import argparse
import glob
import os


def process_run_dir(run_dir):
    pdb = os.path.join(run_dir, "minimized.pdb")
    dcd = os.path.join(run_dir, "prod.dcd")
    out_dcd = os.path.join(run_dir, "prod_centered.dcd")
    # out_pdb = os.path.join(run_dir, "prod_frame0_centered.pdb")

    traj = md.load(dcd, top=pdb)
    protein_idx = traj.topology.select("protein")
    traj.image_molecules(inplace=True)

    # 2) recenter to protein COM per frame
    print("Centering on protein COM...")
    for i, frame in enumerate(traj):
        xyz = frame.xyz[0]
        protein_xyz = xyz[protein_idx]
        com = protein_xyz.mean(axis=0)

        # use box center when available, otherwise center at origin
        box_lengths = frame.unitcell_lengths[0] if frame.unitcell_lengths is not None else None
        has_box = box_lengths is not None and np.any(box_lengths > 0)
        if has_box:
            box_center = box_lengths / 2.0
            shift = box_center - com
        else:
            shift = -com

        frame.xyz[0] = xyz + shift

    traj.image_molecules(inplace=True)
    traj.save(out_dcd)


if __name__ == "__main__":
    for i in range(1,5):
        run_dir = f"../runs_with_equil/ec03_1_runs/run_{i}/run"
        process_run_dir(run_dir)