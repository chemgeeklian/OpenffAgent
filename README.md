```
cd openff_agent
conda install -c -y conda-forge pip 'python>=3.12' 'openff-toolkit-base>=0.17.1' rustworkx rdkit openmm pyxdg gemmi nglview
conda install -y -c conda-forge pdbfixer chemper ambertools # ambertool is needed to add partial charge to PTM
conda install openff-interchange -c conda-forge # required by solvation

pip install -e .
```

order of running scripts for now...

1. core.py
2. struct_fixer.py

build openff force field for PTM (KPI here). Only need to run once for each PTM structure: 

```sh
python process_ptm.py
```

docked ligand anf solvate the protein+ptm_ligand system!

```sh
python dock_ligand.py
python solvate.py
```


# venv on aurora

```sh
cd ../envs
python3 -m venv agent3 --system-site-packages
source /lus/flare/projects/FoundEpidem/xlian/envs/agent3/bin/activate
```

`/lus/flare/projects/FoundEpidem/xlian/envs/agent3/` is the python 3.10 venv to run openmm
If `agent3` is broken, `rm -rf agent3` and `cp -r agent_bkup agent3`
`conda activate agent` activates `/home/xlian/.conda/envs/agent` to run openff...