# Visualize corrected 4mut allosteric effect
# Hsap WT vs 4mut_corrected (D431S+K479E+L495Y+K498T)

reinitialize
cmd.bg_color('white')

# Load
load structures/af3_raw/job1_Hsap_WT/cgas_CT_200-554.pdb, wt
load structures/af3_raw/hsap_cgas_4mut_corrected/model_0.pdb, mut

# Color
as cartoon
color gray80, wt
color salmon, mut

# 4mut sites as blue spheres
show spheres, (wt or mut) and resi 431+479+495+498 and name CA
set sphere_scale, 1.0, (wt or mut) and resi 431+479+495+498 and name CA
color blue, (wt or mut) and resi 431+479+495+498 and name CA

# Displaced N-terminal region as red/green spheres
show spheres, (wt or mut) and resi 211+212+213+214+215+216+217+218+219 and name CA
set sphere_scale, 1.2, (wt or mut) and resi 211+212+213+214+215+216+217+218+219 and name CA
color green, wt and resi 211+212+213+214+215+216+217+218+219 and name CA
color red, mut and resi 211+212+213+214+215+216+217+218+219 and name CA

# Zoom to N-terminal displaced region
zoom resi 200-230

set cartoon_fancy_helices, 1
set antialias, 2
set depth_cue, 0
set orthoscopic, 1

ray 2400, 1800
png figures/hsap_corrected_4mut_allostery.png, dpi=300

# Now Hgal overlay
reinitialize
cmd.bg_color('white')

load structures/af3_raw/job3_Hgal_WT/cgas_CT_200-554.pdb, hgal_wt
load structures/af3_raw/hgal_cgas_4mut_rev_corrected/model_0.pdb, hgal_mut

as cartoon
color gray80, hgal_wt
color salmon, hgal_mut

# 4mut_rev sites
show spheres, (hgal_wt or hgal_mut) and resi 463+511+527+530 and name CA
set sphere_scale, 1.0
color blue, (hgal_wt or hgal_mut) and resi 463+511+527+530 and name CA

# Displaced region
show spheres, (hgal_wt or hgal_mut) and resi 208+209+210+211+212+213 and name CA
set sphere_scale, 1.2
color green, hgal_wt and resi 208+209+210+211+212+213 and name CA
color red, hgal_mut and resi 208+209+210+211+212+213 and name CA

zoom resi 200-260

set cartoon_fancy_helices, 1
set antialias, 2
set depth_cue, 0

ray 2400, 1800
png figures/hgal_corrected_4mut_rev_allostery.png, dpi=300

print("Done")
