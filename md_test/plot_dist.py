#!/usr/bin/env python3
"""
Compute distance from prod.dcd and compare WT vs EC03_1.
"""

from pathlib import Path
import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt

# ============================================================
# CONFIG
# ============================================================
ROOT = Path("/home/xlian/OpenffAgent")
SYSTEMS = {
    "WT": ROOT / "wt_runs",
    "EC03_1": ROOT / "ec03_1_runs",
}

REPEATS = 4
FRAME_DT_PS = 10.0
OUT = ROOT / "plots/WT_vs_EC03_1_prod_distance.png"

# atom selection (CONFIRM NAMES)
SEL_UNK = "resn UNK and name C1x"
SEL_KPI = "resn KPI and name C1"

# ============================================================
# LOAD + COMPUTE
# ============================================================
plt.figure(figsize=(8, 5))

for label, base in SYSTEMS.items():
    color = "tab:blue" if label == "WT" else "tab:red"

    for i in range(1, REPEATS + 1):
        run = base / f"run_{i}" / "run"
        pdb = run / "minimized.pdb"
        dcd = run / "prod.dcd"

        traj = md.load(dcd, top=pdb)

        idx1 = traj.topology.select(SEL_UNK)
        idx2 = traj.topology.select(SEL_KPI)

        if len(idx1) != 1 or len(idx2) != 1:
            raise RuntimeError(f"{label} run {i}: atom selection failed")

        dist = md.compute_distances(traj, [(idx1[0], idx2[0])]).squeeze()
        t = np.arange(len(dist)) * FRAME_DT_PS / 1000.0

        plt.plot(
            t, dist,
            color=color, alpha=0.4,
            label=label if i == 1 else None
        )

# ============================================================
# PLOT
# ============================================================
plt.xlabel("Time (ns)")
plt.ylabel("Distance (Å)")
plt.title("WT vs EC03_1 – PROD distance (from DCD)")
plt.legend()
plt.tight_layout()

OUT.parent.mkdir(exist_ok=True)
plt.savefig(OUT, dpi=300)
plt.close()

print("✅ saved:", OUT)
