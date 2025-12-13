import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


def compute_distance(dcd, top):
    traj = md.load(dcd, top=top)

    unk_idx = traj.topology.select("resn UNK and name C1x")
    kpi_idx = traj.topology.select("resn KPI and name C1")

    if len(unk_idx) != 1 or len(kpi_idx) != 1:
        raise ValueError(f"Found atoms: UNK={unk_idx}, KPI={kpi_idx}")

    i, j = unk_idx[0], kpi_idx[0]
    dist = md.compute_distances(traj, [(i, j)])[:, 0]  # nm
    return dist


systems = {
    "WT": Path("/home/xlian/OpenffAgent/wt_tmp/run"),
    "DgoA*": Path("/home/xlian/OpenffAgent/tmp/run"),
}

results = {}

for name, run_dir in systems.items():
    dcd = run_dir / "prod.dcd"
    top = run_dir / "minimized.pdb"
    out = run_dir / "distance_UNK_C1x_to_KPI_C1.txt"

    print(f"▶ Processing {name}")
    dist = compute_distance(dcd, top)
    np.savetxt(out, dist)
    results[name] = dist


# =========================
# Plot distance vs step
# =========================
plt.figure(figsize=(8, 4))

# Convert frame index to time in nanoseconds.
# Assumption: 1 frame = 10 ps -> 10 ps = 0.01 ns
frame_dt_ps = 10.0
frame_dt_ns = frame_dt_ps / 1000.0

for name, dist in results.items():
    times_ns = np.arange(len(dist)) * frame_dt_ns
    plt.plot(
        times_ns,
        dist,
        label=name,
        linewidth=1.2,
        alpha=0.9,
    )

plt.xlabel("Time (ns)")
plt.ylabel("Distance KPI-C1 LIG-C1X (nm)")
plt.legend()
plt.tight_layout()
plt.savefig("UNK_C1x_to_KPI_C1_distance_vs_time_ns.png", dpi=300)
plt.show()
