#!/usr/bin/env python3
"""
Fix ff19SB CMAP residue-specificity in GROMACS topology.

Problem: parmed converts AMBER prmtop to GROMACS .top, but GROMACS's [cmaptypes]
matching is based on atom types, not residue names. ff19SB has 14 residue-specific
CMAP types, but the converted .top has only 1 (all CA atoms use type 'XC').

Solution:
1. Read AMBER prmtop to get the 14 CMAP grid definitions and the residue->type mapping.
2. For each CA atom in the GROMACS [atoms] section, rename its atom type from 'XC'
   to 'XC{n}' where n is the CMAP type index (0-13) for that residue.
3. In [atomtypes], duplicate the 'XC' entry for each new type 'XC{n}'.
4. In [cmaptypes], replace the single 'C N XC C N' entry with 14 entries
   'C N XC{n} C N', each with the correct grid values from AMBER.
"""
import argparse
import re
from pathlib import Path
from collections import defaultdict

import parmed as pmd


# AMBER CMAP type -> residue name mapping (from analysis of Hsap_WT.prmtop)
CMAP_TYPE_RESIDUES = {
    0: {"GLY"},
    1: {"ALA"},
    2: {"SER"},
    3: {"THR"},
    4: {"ILE", "VAL"},
    5: {"ASN"},
    6: {"ASP"},
    7: {"GLU"},
    8: {"ARG"},
    9: {"LYS"},
    10: {"HID", "HIE", "LEU", "MET", "PHE", "TRP", "TYR"},
    11: {"PRO"},
    12: {"GLN"},
    13: {"CYS"},
}

# Reverse mapping: residue name -> CMAP type index
RES_TO_CMAP_TYPE = {}
for type_idx, residues in CMAP_TYPE_RESIDUES.items():
    for res in residues:
        RES_TO_CMAP_TYPE[res] = type_idx


def get_cmap_data_from_amber(prmtop_path):
    """Extract CMAP grid values and residue mapping from AMBER prmtop."""
    p = pmd.load_file(prmtop_path)

    # Read grid values for each CMAP type
    cmap_grids = {}
    for i, ct in enumerate(p.cmap_types):
        grid_vals = list(ct.grid)  # 576 values (24x24)
        cmap_grids[i] = grid_vals

    # Map each CMAP term to its type index and center residue
    term_mapping = []  # list of (res_name, res_idx, type_idx)
    for cmap in p.cmaps:
        type_idx = p.cmap_types.index(cmap.type)
        ca_atom = p.atoms[cmap.atom3.idx]
        term_mapping.append((
            ca_atom.residue.name,
            ca_atom.residue.idx,
            type_idx,
        ))

    # Verify mapping is consistent
    res_to_type = defaultdict(set)
    for res_name, res_idx, type_idx in term_mapping:
        res_to_type[res_name].add(type_idx)

    # Build residue index -> type_idx mapping (for the first chain)
    residx_to_type = {}
    for res_name, res_idx, type_idx in term_mapping:
        if res_idx not in residx_to_type:
            residx_to_type[res_idx] = type_idx

    return cmap_grids, residx_to_type


def format_cmaptype_line(atom_types, funct, resolution, grid_vals):
    """Format a single cmaptypes entry in GROMACS format."""
    # GROMACS format: at1 at2 at3 at4 at5 funct res1 res2 [grid values...]
    # Grid values are space-separated, can span multiple lines with backslash continuation
    header = f"{atom_types[0]:<6} {atom_types[1]:<6} {atom_types[2]:<6} {atom_types[3]:<6} {atom_types[4]:<6} {funct:>6} {resolution:>6} {resolution:>6}\\"
    lines = [header]
    # Write 576 grid values, 10 per line (common GROMACS convention)
    for i in range(0, len(grid_vals), 10):
        chunk = grid_vals[i:i+10]
        line = " ".join(f"{v:.17g}" for v in chunk)
        if i + 10 < len(grid_vals):
            line += "\\"
        lines.append(line)
    return "\n".join(lines)


def fix_gromacs_top(top_path, prmtop_path, out_top_path):
    """Fix CMAP in GROMACS topology."""
    top_path = Path(top_path)
    out_top_path = Path(out_top_path)

    print(f"Reading AMBER CMAP data from: {prmtop_path}")
    cmap_grids, residx_to_type = get_cmap_data_from_amber(prmtop_path)
    print(f"  Found {len(cmap_grids)} CMAP types")

    print(f"Reading GROMACS topology: {top_path}")
    raw = top_path.read_text()
    
    # Fix compact format: parmed sometimes generates extremely long lines
    # (>4095 chars) which GROMACS cannot read. Split into proper format.
    # First, split by semicolons that separate merged lines (common in parmed output)
    # Then split by newlines normally
    if len(raw) < 10000000:  # sanity check
        # Try to fix the compact format by inserting newlines before section headers
        # and after semicolons that end comments
        import re
        # Split extremely long lines (>4000 chars) at logical boundaries
        fixed_lines = []
        for line in raw.splitlines():
            if len(line) > 4000:
                # This is likely a merged line. Try to split at atom boundaries.
                # GROMACS [atoms] format: nr type resnr residue atom cgnr charge mass
                # Look for patterns like "  num  TYPE  num  RES  ATOM  num  charge  mass"
                # and split between them
                parts = re.split(r'(\s+\d+\s+\S+\s+\d+\s+\S+\s+\S+\s+\d+\s+-?\d+\.\d+\s+\d+\.\d+)', line)
                # Reconstruct properly split lines
                current = ""
                for part in parts:
                    if not part.strip():
                        continue
                    if re.match(r'^\s+\d+\s+\S+\s+\d+\s+\S+\s+\S+\s+\d+\s+-?\d+\.\d+\s+\d+\.\d+', part):
                        if current.strip():
                            fixed_lines.append(current)
                        current = part
                    else:
                        current += part
                if current.strip():
                    fixed_lines.append(current)
            else:
                fixed_lines.append(line)
        lines = fixed_lines
    else:
        lines = raw.splitlines()

    # Parse sections (handle compact format where [section] may share line with data)
    import re
    section_ranges = {}
    current_section = None
    for i, line in enumerate(lines):
        m = re.match(r'^\[\s*(\w+)\s*\]', line)
        if m:
            sec_name = m.group(1)
            if current_section is not None:
                section_ranges[current_section] = (section_ranges[current_section][0], i)
            current_section = sec_name
            section_ranges[current_section] = (i, None)
    # Close last section
    if current_section is not None:
        section_ranges[current_section] = (section_ranges[current_section][0], len(lines))

    # --- Step 1: Fix [atoms] section ---
    # Find CA atoms and rename their type from 'XC' to 'XC{n}'
    atoms_start, atoms_end = section_ranges.get("atoms", (None, None))
    if atoms_start is None:
        raise ValueError("No [atoms] section found in topology")

    # We need to know which moleculetype each [atoms] belongs to
    # GROMACS top has multiple [atoms] sections, one per [moleculetype]
    # We need to track which moleculetype is active

    # Find all [moleculetype] sections and their associated [atoms]
    moltypes = []
    moltype_indices = []
    for i, line in enumerate(lines):
        if line.strip() == "[ moleculetype ]":
            moltype_indices.append(i)

    # For each moleculetype, find its [atoms] section
    atom_sections = []
    for idx in moltype_indices:
        # Find the [atoms] section after this moleculetype
        for j in range(idx + 1, len(lines)):
            if lines[j].strip() == "[ atoms ]":
                atom_sections.append(j)
                break

    print(f"  Found {len(atom_sections)} [atoms] sections (molecules)")

    # Build a mapping from (molecule_idx, atom_index_in_molecule) to cmap_type
    # We need to track residue sequence within each molecule

    new_lines = list(lines)
    modifications = 0

    for mol_idx, atoms_sec_idx in enumerate(atom_sections):
        # Find the end of this [atoms] section
        atoms_end_idx = len(lines)
        for j in range(atoms_sec_idx + 1, len(lines)):
            if lines[j].strip().startswith("[") and lines[j].strip().endswith("]"):
                atoms_end_idx = j
                break

        # Parse atoms in this section
        # Format: nr type resnr residu atom cgnr charge mass [; comment]
        current_resnr = None
        current_resname = None
        res_counter = 0  # 0-based residue index within this molecule

        for line_idx in range(atoms_sec_idx + 1, atoms_end_idx):
            line = new_lines[line_idx]
            stripped = line.strip()
            if not stripped or stripped.startswith(";"):
                continue

            parts = stripped.split()
            if len(parts) < 8:
                continue

            atom_name = parts[4]
            if atom_name == "CA":
                res_name = parts[3]
                res_nr = int(parts[2])

                # Determine CMAP type for this residue
                cmap_type = RES_TO_CMAP_TYPE.get(res_name)
                if cmap_type is None:
                    # Try common AMBER/HIS variants
                    if res_name in ("HIS", "HIP", "HIN"):
                        cmap_type = RES_TO_CMAP_TYPE.get("HIE") or RES_TO_CMAP_TYPE.get("HID")
                    elif res_name == "CYX":
                        cmap_type = RES_TO_CMAP_TYPE.get("CYS")
                    elif res_name == "ASH":
                        cmap_type = RES_TO_CMAP_TYPE.get("ASP")
                    elif res_name == "GLH":
                        cmap_type = RES_TO_CMAP_TYPE.get("GLU")
                    elif res_name == "LYN":
                        cmap_type = RES_TO_CMAP_TYPE.get("LYS")

                if cmap_type is not None:
                    old_type = parts[1]
                    new_type = f"XC{cmap_type}"
                    if old_type == "XC":
                        parts[1] = new_type
                        new_lines[line_idx] = " ".join(parts) + "\n"
                        modifications += 1
                    elif old_type.startswith("XC") and not old_type.startswith("XC0"):
                        # Already partially fixed? Skip
                        pass

    print(f"  Modified {modifications} CA atom types in [atoms]")

    # --- Step 2: Fix [atomtypes] section ---
    # Find 'XC' entry and duplicate it for XC0-XC13
    at_start, at_end = section_ranges.get("atomtypes", (None, None))
    xc_line_idx = None
    xc_line = None
    if at_start is not None:
        for i in range(at_start + 1, at_end):
            line = new_lines[i]
            stripped = line.strip()
            if stripped.startswith("XC ") or stripped.startswith("XC\t"):
                xc_line = line
                xc_line_idx = i
                break

    if xc_line_idx is not None:
        # Parse XC line to get parameters
        parts = xc_line.strip().split(";")[0].split()
        # Format: name mass charge ptype sigma epsilon
        # Insert new entries after XC
        insert_lines = []
        for t in range(14):
            new_name = f"XC{t}"
            new_parts = [new_name] + parts[1:]
            insert_lines.append(" ".join(new_parts) + "\n")
        new_lines = new_lines[:xc_line_idx + 1] + insert_lines + new_lines[xc_line_idx + 1:]
        print(f"  Added 14 atom types (XC0-XC13) to [atomtypes]")
    else:
        print("  WARNING: Could not find 'XC' in [atomtypes]")

    # --- Step 3: Replace [cmaptypes] section ---
    ct_start, ct_end = section_ranges.get("cmaptypes", (None, None))
    if ct_start is not None:
        # Find the old single cmaptype entry and replace with 14 entries
        new_cmaptypes = ["[ cmaptypes ]\n", "\n"]
        for type_idx in range(14):
            grid_vals = cmap_grids[type_idx]
            entry = format_cmaptype_line(
                atom_types=["C", "N", f"XC{type_idx}", "C", "N"],
                funct=1,
                resolution=24,
                grid_vals=grid_vals,
            )
            new_cmaptypes.append(entry + "\n\n")

        new_lines = new_lines[:ct_start] + new_cmaptypes + new_lines[ct_end:]
        print(f"  Replaced [cmaptypes] with 14 residue-specific entries")
    else:
        print("  WARNING: No [cmaptypes] section found")

    # Write output
    out_top_path.write_text("".join(new_lines))
    print(f"\nSaved fixed topology: {out_top_path}")
    return True


def main():
    parser = argparse.ArgumentParser(description="Fix ff19SB CMAP in GROMACS topology")
    parser.add_argument("--top", required=True, help="Input GROMACS .top file")
    parser.add_argument("--prmtop", required=True, help="Reference AMBER prmtop for CMAP data")
    parser.add_argument("--out", required=True, help="Output fixed .top file")
    args = parser.parse_args()

    fix_gromacs_top(args.top, args.prmtop, args.out)


if __name__ == "__main__":
    main()
