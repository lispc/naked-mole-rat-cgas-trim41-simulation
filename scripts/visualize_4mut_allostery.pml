# Visualize Hsap WT vs 4mut allosteric effect on N-terminal region
# Shows that 4mut in C-terminus causes large displacement in N-terminal interface region

reinitialize
cmd.set('ray_trace_mode', 0)
cmd.bg_color('white')

# Load structures
load structures/af3_raw/job1_Hsap_WT/cgas_CT_200-554.pdb, hsap_wt
load structures/af3_raw/Hsap_cGAS_4mut/cgas_CT_200-554.pdb, hsap_4mut

# Color scheme
# WT: green, 4mut: magenta
color green, hsap_wt
color magenta, hsap_4mut

# The 4mut residues (C-terminal) - very similar positions
color blue, hsap_wt and resi 463+511+527+530 and name CA
color blue, hsap_4mut and resi 463+511+527+530 and name CA

# N-terminal interface region showing largest displacement
color red, hsap_wt and resi 211-219 and name CA
color orange, hsap_4mut and resi 211-219 and name CA

# Show as cartoon with transparent surface
as cartoon
show spheres, (name CA) and (resi 211+212+213+214+215+216+217+218+219+463+511+527+530)
set sphere_scale, 0.8

# Labels for key residues
label hsap_wt and resi 213 and name CA, "WT 213"
label hsap_4mut and resi 213 and name CA, "4mut 213"
label hsap_wt and resi 463 and name CA, "4mut site"

# View centered on N-terminal displaced region
zoom hsap_wt and resi 200-230

# Set view for publication
set cartoon_transparency, 0.3
set cartoon_fancy_helices, 1
set cartoon_fancy_sheets, 1
set antialias, 2
set line_smooth, 1
set depth_cue, 0

# Ray trace
ray 2400, 1800
png figures/hsap_4mut_allostery_nterm.png, dpi=300

# Now zoom to C-terminal 4mut region to show they're similar
zoom hsap_wt and resi 460-535
hide labels
label hsap_wt and resi 463 and name CA, "463"
label hsap_wt and resi 511 and name CA, "511"
label hsap_wt and resi 527 and name CA, "527"
label hsap_wt and resi 530 and name CA, "530"

ray 2400, 1800
png figures/hsap_4mut_active_site_similar.png, dpi=300

print("Done: figures/hsap_4mut_allostery_nterm.png")
print("Done: figures/hsap_4mut_active_site_similar.png")
