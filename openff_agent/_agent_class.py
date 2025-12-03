try:
    from openff_agent.run_chai import fold
except:
    print("Could not import fold from openff_agent.run_chai")
from pathlib import Path
from openff_agent import utils
from openff.toolkit import Molecule
import openff_agent.utils as utils


class OpenFFAgent:
    """Main agent for protein structure processing and force field preparation."""
    
    def __init__(self, base_output_dir: Path, ph: float = 7.4):
        """Initialize the agent with configuration."""
        self.base_output_dir = base_output_dir
        self.ph = ph
        self.fasta_jobs = {}
        self.results = {}
    
    # Configuration & Setup
    def register_fasta_job(self, tag: str, fasta_file: Path):
        """Register a FASTA file for processing."""
        self.fasta_jobs[tag] = fasta_file
    
    # Structure Preparation
    def fold_proteins(self):
        """Run CHAI folding for all registered FASTA files."""
        for tag, fasta_file in self.fasta_jobs.items():
            output_dir = self.base_output_dir / tag
            output_dir.mkdir(parents=True, exist_ok=True)
            result = fold(fasta=fasta_file, output_dir=output_dir)
            self.results[tag] = result
    
    def fix_pdb(self, tag: str, resid_to_rm_atom: dict = None):
        """Fix PDB structure (missing atoms, hydrogens, etc.)."""
        pdb_file = self.results[tag]["pdb"]
        fixed_file = pdb_file.replace(".pdb", "_fixed.pdb")
        utils.fix_pdb(pdb_file, fixed_file, resid_to_rm_atom=resid_to_rm_atom, ph=self.ph)
        return fixed_file
    
    # Ligand/KPI Integration
    def load_ligand_smiles(self, smiles_file: Path) -> str:
        """Load ligand SMILES from file."""
        return smiles_file.read_text().strip()
    
    def create_molecule_with_ligand(self, tag: str, fixed_pdb: str, smiles: str):
        """Combine protein and ligand into single molecule."""
        mol = Molecule.from_pdb_and_smiles(fixed_pdb, smiles, allow_undefined_stereo=True)
        self.results[tag]["molecule"] = mol
        return mol
    
    # Run complete pipeline
    def run_pipeline(self, ligand_smiles_file: Path = None):
        """Execute full workflow."""
        self.fold_proteins()
        for tag in self.fasta_jobs.keys():
            self.fix_pdb(tag)
        if ligand_smiles_file:
            smiles = self.load_ligand_smiles(ligand_smiles_file)
            for tag in self.fasta_jobs.keys():
                self.create_molecule_with_ligand(tag)