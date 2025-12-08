import MDAnalysis as mda

#u = mda.Universe("../EC03_1_sim/system_rm_CL.pdb")
#u = mda.Universe("dgoasystem.pdb")
u = mda.Universe("/lus/flare/projects/FoundEpidem/xlian/OpenffAgent/md_test/run/system.pdb")

# Use built-in guesses for residue charges (neutral termini assumed)
total_charge = 0.0
for res in u.residues:
    resname = res.resname.strip().upper()
    if resname in ["ASP", "GLU"]: total_charge -= 1
    elif resname in ["LYS", "ARG"]: total_charge += 1
    # optional: add metal ions
    elif resname in ["NA"]: total_charge += 1
    elif resname in ["CL"]: total_charge -= 1

total_charge -= 3 # pyruvate -1, E4P -2

print("Estimated total charge:", total_charge)