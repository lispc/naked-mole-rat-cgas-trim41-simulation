# PyMOL script: Visualize Hgal best_pose interface
# Shows that the actual interface is at cGAS N-terminus (resi 213-247),
# while the "active site" residues (463/511/527/530) are 30-39A away.
#
# Run: pymol -cq scripts/visualize_hgal_interface.pml
# Or:  pymol scripts/visualize_hgal_interface.pml

reinitialize
bg_color white
set depth_cue, 0
set orthoscopic, 1
set cartoon_fancy_helices, 1
set cartoon_smooth_loops, 1

# Load processed PDB (chains already split: A=TRIM41, B=cGAS)
load data/md_runs/Hgal_domain/Hgal_domain_processed.pdb, complex

# Color scheme
# TRIM41: light gray
color gray70, complex and chain A

# cGAS N-terminal region (200-300): light blue
color marine, complex and chain B and resi 200-300

# cGAS C-terminal / active site region (400-554): pale cyan
color palecyan, complex and chain B and resi 400-554

# cGAS middle (301-399): very pale blue
color lightblue, complex and chain B and resi 301-399

# Show cartoon
show cartoon, complex

# Actual interface residues (from 5A cutoff analysis)
select interface_cgas, complex and chain B and resi 213+242+247
select interface_trim41, complex and chain A and resi 571+573+598
select interface_all, interface_cgas or interface_trim41

# Show interface as sticks
show sticks, interface_all
set stick_radius, 0.3, interface_all
color red, interface_cgas
color orange, interface_trim41

# Active site residues (the 4 "mutations")
select active_site, complex and chain B and resi 463+511+527+530 and name CA
show spheres, active_site
set sphere_scale, 1.5, active_site
color blue, active_site

# Label active site residues
label active_site and resi 463, "S463"
label active_site and resi 511, "E511"
label active_site and resi 527, "Y527"
label active_site and resi 530, "T530"
set label_color, blue, active_site
set label_size, 20

# Label interface residues
label interface_cgas and resi 213 and name CA, "I213"
label interface_cgas and resi 242 and name CA, "R242"
label interface_cgas and resi 247 and name CA, "E247"
set label_color, red, interface_cgas

# Measure distances from active sites to nearest TRIM41 atom
# These will be shown as dashed lines
# (PyMOL may not have exact atoms; we'll use CA-CA approximations)

# Zoom to show whole complex
zoom complex, 20

# Ray trace and save
set ray_trace_mode, 0
set opaque_background, 1
set ray_shadows, 0

# View 1: Full complex overview
view 0, store
ray 2400, 1800
png figures/hgal_interface_overview.png, dpi=300

# View 2: Zoom to interface region
zoom interface_all, 10
ray 2400, 1800
png figures/hgal_interface_closeup.png, dpi=300

# View 3: Zoom to active site region (showing they are NOT at interface)
zoom active_site, 15
ray 2400, 1800
png figures/hgal_interface_active_site.png, dpi=300

print("Figures saved to figures/hgal_interface_*.png")
print("Key finding: Interface at N-terminus (213/242/247), NOT at active site (463/511/527/530)")
print("Active site residues are 30-39A away from TRIM41")

quit
