#!/usr/bin/env python3
"""
Compute average RMSD from prod.dcd and compare WT vs EC03_1.
(mean ± std over repeats)
"""

from pathlib import Path
import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt

# ============================================================
# CONFIG
# ============================================================
ROOT = Path("/home/xlian/OpenffAgent/runs_no_equil")
SYSTEMS = {
    "WT": ROOT / "wt_runs",
    "EC03_1": ROOT / "ec03_1_runs",
}

REPEATS = 4
FRAME_DT_PS = 10.0
OUT = ROOT / "plots/WT_vs_EC03_1_prod_rmsd_avg.png"

# atom selection for RMSD calculation
SEL_ATOMS = "protein and name CA"  # C-alpha atoms for protein backbone RMSD

# ============================================================
# LOAD + COMPUTE
# ============================================================
plt.figure(figsize=(8, 5))

for label, base in SYSTEMS.items():
    color = "tab:blue" if label == "WT" else "tab:red"

    all_rmsd = []

    for i in range(1, REPEATS + 1):
        run = base / f"run_{i}" / "run"
        pdb = run / "minimized.pdb"
        dcd = run / "prod.dcd"

        traj = md.load(dcd, top=pdb)

        # Select atoms for RMSD calculation
        atom_indices = traj.topology.select(SEL_ATOMS)

        if len(atom_indices) == 0:
            raise RuntimeError(f"{label} run {i}: atom selection failed")

        # Compute RMSD relative to first frame
        rmsd = md.rmsd(traj, traj, 0, atom_indices=atom_indices)

        all_rmsd.append(rmsd)

    # --------------------------------------------------------
    # Stack & average
    # --------------------------------------------------------

    min_len = min(len(r) for r in all_rmsd)
    all_rmsd = np.array([r[:min_len] for r in all_rmsd])

    mean_rmsd = all_rmsd.mean(axis=0)
    std_rmsd = all_rmsd.std(axis=0)

    t = np.arange(mean_rmsd.shape[0]) * FRAME_DT_PS / 1000.0

    # --------------------------------------------------------
    # Plot
    # --------------------------------------------------------
    plt.plot(t, mean_rmsd, color=color, label=label, linewidth=2)
    plt.fill_between(
        t,
        mean_rmsd - std_rmsd,
        mean_rmsd + std_rmsd,
        color=color,
        alpha=0.25
    )

# ============================================================
# PLOT
# ============================================================
plt.xlabel("Time (ns)")
plt.ylabel("RMSD (nm)")
plt.title("WT vs EC03_1 – PROD RMSD (mean ± std)")
plt.legend()
plt.tight_layout()

OUT.parent.mkdir(exist_ok=True)
plt.savefig(OUT, dpi=300)
plt.close()

print("Saved:", OUT)

