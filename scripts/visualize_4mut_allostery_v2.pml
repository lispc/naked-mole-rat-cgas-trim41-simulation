# Hsap WT vs 4mut allosteric effect visualization v2
# Clean overlay showing N-terminal displacement

reinitialize
cmd.set('ray_trace_mode', 0)
cmd.bg_color('white')

# Load
load structures/af3_raw/job1_Hsap_WT/cgas_CT_200-554.pdb, wt
load structures/af3_raw/Hsap_cGAS_4mut/cgas_CT_200-554.pdb, mut

# Hide everything initially
hide everything

# Show as cartoon
as cartoon

# Color by chain
# WT: gray, 4mut: salmon (overlay)
color gray80, wt
color salmon, mut

# Highlight displaced N-terminal region (200-230) as thicker cartoon
set cartoon_transparency, 0.0, wt and resi 200-230
set cartoon_transparency, 0.0, mut and resi 200-230
set cartoon_transparency, 0.7, wt and resi 300-554
set cartoon_transparency, 0.7, mut and resi 300-554

# Highlight specific displaced residues as spheres
show spheres, (wt or mut) and resi 211+212+213+214+215+216+217+218+219 and name CA
set sphere_scale, 1.2, (wt or mut) and resi 211+212+213+214+215+216+217+218+219 and name CA

# Color spheres: WT green, 4mut red
color green, wt and resi 211+212+213+214+215+216+217+218+219 and name CA
color red, mut and resi 211+212+213+214+215+216+217+218+219 and name CA

# 4mut sites as blue spheres
show spheres, (wt or mut) and resi 463+511+527+530 and name CA
set sphere_scale, 1.0, (wt or mut) and resi 463+511+527+530 and name CA
color blue, (wt or mut) and resi 463+511+527+530 and name CA

# Zoom to N-terminal displaced region
zoom resi 200-230

# Publication settings
set cartoon_fancy_helices, 1
set antialias, 2
set depth_cue, 0
set orthoscopic, 1

ray 2400, 1800
png figures/hsap_4mut_allostery_clean.png, dpi=300

# Side-by-side view
reinitialize
cmd.bg_color('white')

# WT
load structures/af3_raw/job1_Hsap_WT/cgas_CT_200-554.pdb, wt
as cartoon
color gray80
color green, resi 211-219
color blue, resi 463+511+527+530
show spheres, resi 211+212+213+214+215+216+217+218+219+463+511+527+530 and name CA
set sphere_scale, 1.0
zoom resi 200-230
set_view (\
    -0.67, 0.48, -0.56,\
    -0.73, -0.31, 0.61,\
    0.12, 0.82, 0.56,\
    0.00, 0.00, -80.00,\
    15.00, 5.00, -5.00,\
    60.00, 100.00, -20.00)

ray 2400, 1800
png figures/hsap_wt_nterm.png, dpi=300

# 4mut
reinitialize
cmd.bg_color('white')
load structures/af3_raw/Hsap_cGAS_4mut/cgas_CT_200-554.pdb, mut
as cartoon
color gray80
color red, resi 211-219
color blue, resi 463+511+527+530
show spheres, resi 211+212+213+214+215+216+217+218+219+463+511+527+530 and name CA
set sphere_scale, 1.0
zoom resi 200-230
set_view (\
    -0.67, 0.48, -0.56,\
    -0.73, -0.31, 0.61,\
    0.12, 0.82, 0.56,\
    0.00, 0.00, -80.00,\
    15.00, 5.00, -5.00,\
    60.00, 100.00, -20.00)

ray 2400, 1800
png figures/hsap_4mut_nterm.png, dpi=300

print("Done")
