"""PSF (Protein Structure File) parser — produces bonds, angles, dihedrals, impropers."""


def from_psf(psf_path: str) -> dict:
    """Parse a PSF file into a system.json dict (force_field section only)."""
    bonds = []
    angles = []
    dihedrals = []
    impropers = []
    section = None
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

            if section == "bonds":
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
    return {"force_field": force_field}
