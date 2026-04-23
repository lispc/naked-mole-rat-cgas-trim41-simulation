# PyMOL script: Compare Rosetta docking with LightDock reference
# Run with: pymol compare_visualization.pml

# Load reference (LightDock)
load input.pdb, lightdock_ref
hide everything, lightdock_ref
show cartoon, lightdock_ref
color cyan, lightdock_ref and chain A
color marine, lightdock_ref and chain B

# Load best Rosetta decoy
load output_global/input_0003.pdb, rosetta_best
hide everything, rosetta_best
show cartoon, rosetta_best
color orange, rosetta_best and chain A
color red, rosetta_best and chain B

# Align Rosetta to reference (chain B = cGAS CTD as ligand)
align rosetta_best and chain B, lightdock_ref and chain B

# Show interface residues (within 5Å)
select interface_ref, (lightdock_ref and chain A) within 5 of (lightdock_ref and chain B)
select interface_ros, (rosetta_best and chain A) within 5 of (rosetta_best and chain B)
show sticks, interface_ref
show sticks, interface_ros
color yellow, interface_ref
color magenta, interface_ros

# Labels
title "Rosetta (orange/red) vs LightDock (cyan/marine)\nBest decoy: input_0003, RMSD=2.10Å, I_sc=-23.02 REU"

# Ray trace settings
set ray_trace_mode, 1
set antialias, 2
set cartoon_fancy_helices, 1

print "Visualization loaded. Rosetta best decoy aligned to LightDock reference."
print "Use 'ray 2400,1600' to render high-res image."
