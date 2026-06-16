import numpy as np


def from_pdb(pdb_path: str) -> dict:
    """Parse a PDB file into a system.json dict."""
    atoms = []
    positions = []
    box_size = 3.0
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("CRYST1"):
                parts = line.split()
                box_size = float(parts[1]) if len(parts) > 1 else 3.0
            elif line.startswith("ATOM") or line.startswith("HETATM"):
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                positions.append([x, y, z])
                elem = line[76:78].strip() or line[12:16].strip()
                atoms.append({
                    "mass": _guess_mass(elem),
                    "charge": 0.0,
                    "type": elem,
                })
    return {
        "cell": {"positions": positions, "box_size": box_size},
        "atoms": {k: [a[k] for a in atoms] for k in atoms[0]},
        "force_field": {},
    }


def _guess_mass(elem: str) -> float:
    masses = {
        "H": 1.008, "HE": 4.003, "C": 12.011, "N": 14.007,
        "O": 15.999, "F": 18.998, "NA": 22.990, "MG": 24.305,
        "P": 30.974, "S": 32.065, "CL": 35.450, "K": 39.098,
        "CA": 40.078, "AR": 39.948,
    }
    return masses.get(elem.upper().strip(), 12.011)
