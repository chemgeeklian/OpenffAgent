from openff_agent.struct_fixer import prepare_structure
from pathlib import Path

single_seq_fasta_file = '../data/chai_input.fasta'
output_dir = Path('../tmp')
chai_output = output_dir / "pred.model_idx_0.cif"

chai_output_pdb = chai_output.with_suffix(".pdb")

if False:
    from openff_agent.run_chai import fold
    result = fold(
        fasta=Path(single_seq_fasta_file),
        output_dir=Path(output_dir),
        constraint=None,
        fold_algo="chai-esm"
    )

    from openff_agent.utils import cif2pdb
    cif2pdb(str(chai_output))

prepare_structure(chai_output_pdb, chai_output.with_name(chai_output.stem + "_H.pdb"), ph=7.4)