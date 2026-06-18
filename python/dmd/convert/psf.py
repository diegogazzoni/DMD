from __future__ import annotations
"""PSF (Protein Structure File) parser — produces bonds, angles, dihedrals, impropers + atom types + exclusions."""


def _exclusion_pairs(
    bonds: list[dict],
    angles: list[dict] | None = None,
    dihedrals: list[dict] | None = None,
    scale_14: float = 1.0,
) -> list[list[int]]:
    """Build exclusion pair list from topology data.

    1-2 pairs (bonds) and 1-3 pairs (angles) are fully excluded.
    1-4 pairs (dihedrals) are fully excluded if scale_14=1.0.
    """
    excl = set()
    for b in bonds:
        i, j = b["i"], b["j"]
        excl.add((i, j) if i < j else (j, i))
    if angles:
        for a in angles:
            i, k = a["i"], a["k"]
            excl.add((i, k) if i < k else (k, i))
    if dihedrals and abs(scale_14 - 1.0) > 1e-9:
        for d in dihedrals:
            i, l = d["i"], d["l"]
            excl.add((i, l) if i < l else (l, i))
    return [list(p) for p in sorted(excl)]


def from_psf(psf_path: str, scale_14: float = 1.0) -> dict:
    """Parse a PSF file into a system.json dict (force_field + atom_types section).

    Returns a dict with:
        force_field: bonds, angles, dihedrals, impropers (indices 0-based)
        atom_types: list of type names per atom (index position = atom index)
        exclusions: list of [i,j] pairs excluded from non-bonded

    Args:
        psf_path: path to the .psf file
        scale_14: 1-4 scaling factor (1.0 = full exclusion, 0.0 = full interaction)
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
                for idx in range(0, len(parts), 2):
                    if idx + 1 < len(parts):
                        bonds.append({
                            "i": int(parts[idx]) - 1,
                            "j": int(parts[idx + 1]) - 1,
                        })
            elif section == "angles":
                parts = stripped.split()
                for idx in range(0, len(parts), 3):
                    if idx + 2 < len(parts):
                        angles.append({
                            "i": int(parts[idx]) - 1,
                            "j": int(parts[idx + 1]) - 1,
                            "k": int(parts[idx + 2]) - 1,
                        })
            elif section == "dihedrals":
                parts = stripped.split()
                for idx in range(0, len(parts), 4):
                    if idx + 3 < len(parts):
                        dihedrals.append({
                            "i": int(parts[idx]) - 1,
                            "j": int(parts[idx + 1]) - 1,
                            "k": int(parts[idx + 2]) - 1,
                            "l": int(parts[idx + 3]) - 1,
                        })
            elif section == "impropers":
                parts = stripped.split()
                for idx in range(0, len(parts), 4):
                    if idx + 3 < len(parts):
                        impropers.append({
                            "i": int(parts[idx]) - 1,
                            "j": int(parts[idx + 1]) - 1,
                            "k": int(parts[idx + 2]) - 1,
                            "l": int(parts[idx + 3]) - 1,
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

    exclusions = _exclusion_pairs(bonds, angles, dihedrals, scale_14)

    result = {"force_field": force_field}
    if atom_types_list:
        result["atom_types"] = atom_types_list
    if exclusions:
        result["exclusions"] = exclusions
    return result
