import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


def compute_distance(dcd, top):
    """Load a trajectory and compute distance between UNK:C1x and KPI:C1.

    Args:
        dcd: path-like to the trajectory file (DCD)
        top: path-like to the topology file (PDB)

    Returns:
        1D numpy array of distances (nm) for each frame.
    """
    traj = md.load(dcd, top=top)

    unk_idx = traj.topology.select("resn UNK and name C1x")
    kpi_idx = traj.topology.select("resn KPI and name C1")

    if len(unk_idx) != 1 or len(kpi_idx) != 1:
        raise ValueError(f"Found atoms: UNK={unk_idx}, KPI={kpi_idx}")

    i, j = unk_idx[0], kpi_idx[0]
    dist = md.compute_distances(traj, [(i, j)])[:, 0]  # nm
    return dist


def process_systems(systems):
    """Process multiple systems: compute distances and save per-system outputs.

    Args:
        systems: mapping of name -> Path (or path-like) to the run directory.

    Returns:
        dict mapping system name to numpy distance arrays.
    """
    results = {}
    for name, run_dir in systems.items():
        run_dir = Path(run_dir)
        dcd = run_dir / "prod.dcd"
        top = run_dir / "minimized.pdb"
        out = run_dir / "distance_UNK_C1x_to_KPI_C1.txt"

        print(f"▶ Processing {name}")
        dist = compute_distance(dcd, top)
        np.savetxt(out, dist)
        results[name] = dist

    return results


def plot_results(results, 
                 frame_dt_ps=10.0, 
                 filename="LIG_C1x_to_KPI_C1_distance_vs_time_ns.png", 
                 show=True):
    """Plot distances vs time for provided results.

    Args:
        results: mapping name -> 1D numpy distances
        frame_dt_ps: time between frames in picoseconds (default 10 ps)
        filename: output filename for saved plot
        show: whether to call `plt.show()` after saving
    """
    plt.figure(figsize=(8, 4))

    # Convert frame index to time in nanoseconds.
    frame_dt_ns = frame_dt_ps / 1000.0

    for name, dist in results.items():
        times_ns = np.arange(len(dist)) * frame_dt_ns
        plt.plot(times_ns, dist, label=name, linewidth=1.2, alpha=0.9)

    plt.xlabel("Time (ns)")
    plt.ylabel("Distance KPI-C1 LIG-C1X (nm)")
    plt.legend()
    plt.tight_layout()
    # Ensure output directory exists
    try:
        fpath = Path(filename)
        if not fpath.parent.exists():
            fpath.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(str(fpath), dpi=300)
    except Exception:
        # Fallback: attempt to save to working directory
        plt.savefig(filename, dpi=300)
    if show:
        plt.show()
    plt.close()


def main():
    systems = {
        "WT": Path("/home/xlian/OpenffAgent/wt_tmp/run"),
        "DgoA*": Path("/home/xlian/OpenffAgent/tmp/run"),
    }

    results = process_systems(systems)
    plot_results(results)


if __name__ == "__main__":
    main()
