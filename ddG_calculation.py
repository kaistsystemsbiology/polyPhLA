from pyrosetta import *
from pyrosetta.toolbox import *
from pyrosetta.rosetta.protocols.relax import FastRelax

### Note: Pyrosequence residue numbers may differ from PDB numbers because they require consecutive numbers for calculation. 

# PyRosetta initialization
init(extra_options="-beta_nov16_cart -in:file:s 1kjg.clean.pdb -use_input_sc -constrain_relax_to_start_coords -ignore_unrecognized_res -relax:coord_constrain_sidechains -relax:ramp_constraints false -relax:cartesian -relax:min_type lbfgs_armijo_nonmonotone")

res_num = 80
origin_res = 'A'
mut_res = 'G'

# PDB loading and energy calculation for the original structure  
pose = pose_from_pdb("../data/phac1437_WT.pdb")
mutate_residue(pose, res_num, origin_res)
original_scorefxn = get_fa_scorefxn()  
original_energy = original_scorefxn(pose)  

# Introduce a single mutation
mutate_residue(pose, res_num, mut_res)

# Optimized structure storage
pose.dump_pdb("phac1437_WT_opt.pdb")
mutated_energy = original_scorefxn(pose)

# ΔΔG calculation
ddG = mutated_energy - original_energy
print(f"The ΔΔG of mutating residue {res_num} to {mut_res} is {ddG:.3f} kcal/mol")