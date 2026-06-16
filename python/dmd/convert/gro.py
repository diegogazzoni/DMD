"""GRO (GROMACS) file parser."""


def from_gro(gro_path: str) -> dict:
    """Parse a GRO file into a system.json dict."""
    with open(gro_path) as f:
        lines = f.readlines()

    n_atoms = int(lines[1].strip())
    atoms = []
    positions = []
    box_line = lines[-1].strip()
    box_parts = box_line.split()
    box_size = float(box_parts[0]) if box_parts else 3.0

    for line in lines[2:2 + n_atoms]:
        resname = line[5:10].strip()
        atomname = line[10:15].strip()
        x = float(line[20:28].strip())
        y = float(line[28:36].strip())
        z = float(line[36:44].strip())
        positions.append([x, y, z])
        atoms.append({
            "mass": _guess_mass(atomname),
            "charge": 0.0,
            "type": atomname,
            "residue": resname,
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
    return masses.get(elem.upper().strip()[:2], 12.011)
