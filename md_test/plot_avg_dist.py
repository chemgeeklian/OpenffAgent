#!/usr/bin/env python3
"""
Compute average distance from prod.dcd and compare WT vs EC03_1.
(mean ± std over repeats)
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
OUT = ROOT / "plots/WT_vs_EC03_1_prod_distance_avg.png"

# atom selection
SEL_UNK = "resn UNK and name C1x"
SEL_KPI = "resn KPI and name C1"

# ============================================================
# LOAD + COMPUTE
# ============================================================
plt.figure(figsize=(8, 5))

for label, base in SYSTEMS.items():
    color = "tab:blue" if label == "WT" else "tab:red"

    all_dist = []

    for i in range(1, REPEATS + 1):
        run = base / f"run_{i}" / "run"
        pdb = run / "minimized.pdb"
        dcd = run / "prod.dcd"

        traj = md.load(dcd, top=pdb)

        idx1 = traj.topology.select(SEL_UNK)
        idx2 = traj.topology.select(SEL_KPI)

        if len(idx1) != 1 or len(idx2) != 1:
            raise RuntimeError(f"{label} run {i}: atom selection failed")

        dist = md.compute_distances(
            traj, [(idx1[0], idx2[0])]
        ).squeeze()

        all_dist.append(dist)

    # --------------------------------------------------------
    # Stack & average
    # --------------------------------------------------------

    min_len = min(len(d) for d in all_dist)
    all_dist = np.array([d[:min_len] for d in all_dist])

    #all_dist = np.array(all_dist)      # (n_repeats, n_frames)

    mean_dist = all_dist.mean(axis=0)
    std_dist = all_dist.std(axis=0)

    t = np.arange(mean_dist.shape[0]) * FRAME_DT_PS / 1000.0

    # --------------------------------------------------------
    # Plot
    # --------------------------------------------------------
    plt.plot(t, mean_dist, color=color, label=label, linewidth=2)
    plt.fill_between(
        t,
        mean_dist - std_dist,
        mean_dist + std_dist,
        color=color,
        alpha=0.25
    )

# ============================================================
# PLOT
# ============================================================
plt.xlabel("Time (ns)")
plt.ylabel("Distance (Å)")
plt.title("WT vs EC03_1 – PROD distance (mean ± std)")
plt.legend()
plt.tight_layout()

OUT.parent.mkdir(exist_ok=True)
plt.savefig(OUT, dpi=300)
plt.close()

print("Saved:", OUT)
