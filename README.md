```
cd openff_agent
conda install -c -y conda-forge pip 'python>=3.12' 'openff-toolkit-base>=0.17.1' rustworkx rdkit openmm pyxdg gemmi
conda install -y -c conda-forge pdbfixer chemper ambertools

pip install -e .
```

order of running scripts for now...

1. core.py
2. struct_fixer.py

build openff force field for PTM (KPI here). Only need to run once for each PTM structure: 

```
python process_ptm.py
```