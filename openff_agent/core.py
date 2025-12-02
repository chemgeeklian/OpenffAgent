from pathlib import Path

# input files
seq_noptm_fasta_file = Path("../data/chai_fasta_noptm.fasta")
seq_ptm_fasta_file   = Path("../data/chai_fasta_ptm.fasta")

# base output folder
base_output = Path("../tmp")

# we want to process both
fasta_jobs = {
    "noptm": seq_noptm_fasta_file,
    "ptm":   seq_ptm_fasta_file,
}

for tag, fasta_file in fasta_jobs.items():

    # each FASTA has its own output subfolder
    output_dir = base_output / tag
    output_dir.mkdir(parents=True, exist_ok=True)

    from openff_agent.run_chai import fold
    print(f"\n🧬 Running CHAI for: {tag} ({fasta_file})")
    result = fold(
        fasta=fasta_file,
        output_dir=output_dir,
        constraint=None,
        fold_algo="chai-esm",
    )

    #from openff_agent.utils import cif2pdb
    #cif2pdb(str(chai_output))

#prepare_structure(chai_output_pdb, chai_output.with_name(chai_output.stem + "_H.pdb"), ph=7.4)