#!/usr/bin/env python
import os
import argparse
from pathlib import Path
from typing import Union

from chai_lab.chai1 import run_inference


def fold(
    fasta: Path,
    output_dir: Path,
    constraint: Union[None, Path] = None,
    fold_algo: str = "chai-esm",
    num_diffn_samples: int = 1,
    device: str = "cuda",
):
    """
    Run CHAI or CHAI-ESM folding inference.
    """
    esm = "esm" in fold_algo.lower()  # Only enable ESM if 'esm' is in fold_algo

    # Ensure correct Path types
    fasta = Path(fasta)
    output_dir = Path(output_dir)
    if constraint is not None:
        constraint = Path(constraint)

    # Auto-create output dir
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"📄 FASTA: {fasta}")
    print(f"📁 Output dir: {output_dir}")
    print(f"🧩 Constraint: {constraint}")
    print(f"🧠 Fold algo: {fold_algo} (use_esm={esm})")
    print(f"🔢 diffN samples: {num_diffn_samples}")
    print(f"🖥️ device: {device}")

    candidates = run_inference(
        fasta_file=fasta,
        output_dir=output_dir,
        num_trunk_recycles=3,
        num_diffn_samples=num_diffn_samples,
        num_diffn_timesteps=200,
        device=device,
        use_esm_embeddings=esm,
        use_msa_server=not esm,
        constraint_path=constraint,
    )

    return candidates


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run CHAI inference with optional custom inputs."
    )

    parser.add_argument(
        "--fasta",
        type=str,
        default="/eagle/projects/FoundEpidem/xlian/IDEAL/nibrin/nibrin.fasta",
        help="Path to input FASTA file.",
    )

    parser.add_argument(
        "--constraint",
        type=str,
        default=None,
        help="Optional constraint file (PDB/JSON).",
    )

    parser.add_argument(
        "--outdir",
        type=str,
        default="/home/xlian/nbn_chai_15samples/",
        help="Output directory.",
    )

    parser.add_argument(
        "--algo",
        type=str,
        default="chai-esm",
        choices=["chai", "chai-esm"],
        help="Use CHAI (MSA) or CHAI-ESM (no MSA).",
    )

    parser.add_argument(
        "--num-diffn-samples",
        type=int,
        default=1,
        help="Number of diffN samples to draw (int).",
    )

    parser.add_argument(
        "--device",
        type=str,
        default="cuda",
        help="Device to run inference on (e.g., 'cpu' or 'cuda').",
    )

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    fold(
        fasta=Path(args.fasta),
        output_dir=Path(args.outdir),
        constraint=Path(args.constraint) if args.constraint else None,
        fold_algo=args.algo,
        num_diffn_samples=args.num_diffn_samples,
        device=args.device,
    )
