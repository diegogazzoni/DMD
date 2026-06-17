"""PSF (Protein Structure File) parser — produces bonds, angles, dihedrals, impropers + atom types."""


def from_psf(psf_path: str) -> dict:
    """Parse a PSF file into a system.json dict (force_field + atom_types section).

    Returns a dict with:
        force_field: bonds, angles, dihedrals, impropers (indices 0-based)
        atom_types: list of type names per atom (index position = atom index)
    """
    bonds = []
    angles = []
    dihedrals = []
    impropers = []
    atom_types_list = []
    section = None
    n_atoms = 0
    atom_line_idx = 0

    with open(psf_path) as f:
        for line in f:
            stripped = line.strip()
            if "!NATOM" in stripped:
                section = "atoms"
                n_atoms = int(stripped.split()[0])
                continue
            elif "!NBOND" in stripped:
                section = "bonds"
                continue
            elif "!NTHETA" in stripped:
                section = "angles"
                continue
            elif "!NPHI" in stripped:
                section = "dihedrals"
                continue
            elif "!NIMPHI" in stripped:
                section = "impropers"
                continue
            elif "!" in stripped and section is not None:
                section = None
                continue

            if section == "atoms":
                if atom_line_idx >= n_atoms:
                    continue
                cols = stripped.split()
                if len(cols) >= 7:
                    atom_type = cols[5]
                elif len(cols) >= 4:
                    atom_type = cols[3]
                else:
                    atom_type = f"UNK{atom_line_idx}"
                atom_types_list.append(atom_type)
                atom_line_idx += 1

            elif section == "bonds":
                parts = stripped.split()
                for i in range(0, len(parts), 2):
                    if i + 1 < len(parts):
                        bonds.append({
                            "i": int(parts[i]) - 1,
                            "j": int(parts[i + 1]) - 1,
                        })
            elif section == "angles":
                parts = stripped.split()
                for i in range(0, len(parts), 3):
                    if i + 2 < len(parts):
                        angles.append({
                            "i": int(parts[i]) - 1,
                            "j": int(parts[i + 1]) - 1,
                            "k": int(parts[i + 2]) - 1,
                        })
            elif section == "dihedrals":
                parts = stripped.split()
                for i in range(0, len(parts), 4):
                    if i + 3 < len(parts):
                        dihedrals.append({
                            "i": int(parts[i]) - 1,
                            "j": int(parts[i + 1]) - 1,
                            "k": int(parts[i + 2]) - 1,
                            "l": int(parts[i + 3]) - 1,
                        })
            elif section == "impropers":
                parts = stripped.split()
                for i in range(0, len(parts), 4):
                    if i + 3 < len(parts):
                        impropers.append({
                            "i": int(parts[i]) - 1,
                            "j": int(parts[i + 1]) - 1,
                            "k": int(parts[i + 2]) - 1,
                            "l": int(parts[i + 3]) - 1,
                        })

    force_field = {}
    if bonds:
        force_field["bonds"] = bonds
    if angles:
        force_field["angles"] = angles
    if dihedrals:
        force_field["dihedrals"] = dihedrals
    if impropers:
        force_field["impropers"] = impropers

    result = {"force_field": force_field}
    if atom_types_list:
        result["atom_types"] = atom_types_list
    return result
