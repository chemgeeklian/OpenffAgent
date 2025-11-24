from run_chai import fold
from pathlib import Path

result = fold(
    fasta=Path("../data/chai_input.fasta"),
    output_dir=Path("../tmp"),
    constraint=None,
    fold_algo="chai-esm"
)