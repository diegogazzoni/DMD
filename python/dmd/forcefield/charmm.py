"""CHARMM PAR/PRM force field parameter file parser."""

import re
from typing import Dict
from dmd.forcefield.base import FFParams
from dmd.forcefield.registry import register

TORAD = 0.017453292519943295  # π/180
KJPERKCAL = 4.184


class CharmmParser:
    """Parse CHARMM PAR/PRM force field parameter files.

    Args:
        convert_units: If True, convert CHARMM units (kcal/mol, Å, degrees)
                       to DMD units (kJ/mol, nm, radians). If False, store
                       values as-is from the PRM file (DMD-friendly units).
                       Default: False (backward compatible with tip3p.prm).
    """

    def __init__(self, convert_units: bool = False):
        self.convert_units = convert_units

    def parse(self, path: str) -> FFParams:
        params = FFParams()
        section = None
        with open(path) as f:
            for line in f:
                stripped = line.strip()
                if not stripped or stripped.startswith("*") or stripped.startswith("!"):
                    continue
                upper = stripped.upper()
                if upper in ("BONDS", "ANGLES", "DIHEDRALS", "IMPROPER", "NONBONDED", "HBOND", "NBFIX", "CMAP"):
                    section = upper
                    continue
                if section == "BONDS":
                    cols = stripped.split()
                    if len(cols) >= 4 and not cols[0].startswith("!"):
                        t1, t2, kb, b0 = cols[0], cols[1], cols[2], cols[3]
                        k = float(kb)
                        r0 = float(b0)
                        if self.convert_units:
                            k *= KJPERKCAL / 0.01  # kcal/mol/Å² → kJ/(mol·nm²)
                            r0 *= 0.1              # Å → nm
                        params.bonds[(t1, t2)] = params.bonds[(t2, t1)] = {
                            "k": k,
                            "r0": r0,
                        }
                elif section == "ANGLES":
                    cols = stripped.split()
                    if len(cols) >= 5 and not cols[0].startswith("!"):
                        t1, t2, t3 = cols[0], cols[1], cols[2]
                        ktheta, theta0 = cols[3], cols[4]
                        key = (t1, t2, t3)
                        kth = float(ktheta)
                        th0 = float(theta0)
                        if self.convert_units:
                            kth *= KJPERKCAL  # kcal/(mol·rad²) → kJ/(mol·rad²)
                            # theta0 also converted below
                        params.angles[key] = {
                            "k_theta": kth,
                            "theta0": th0 * TORAD,  # always degrees → radians
                        }
                elif section == "DIHEDRALS":
                    cols = stripped.split()
                    if len(cols) >= 7 and not cols[0].startswith("!"):
                        t1, t2, t3, t4 = cols[0], cols[1], cols[2], cols[3]
                        kchi, periodicity, delta = cols[4], cols[5], cols[6]
                        key = (t1, t2, t3, t4)
                        kp = float(kchi)
                        ph0 = float(delta)
                        if self.convert_units:
                            kp *= KJPERKCAL  # kcal/mol → kJ/mol
                            ph0 *= TORAD     # degrees → radians
                        entry = {
                            "k_phi": kp,
                            "periodicity": int(periodicity),
                            "phi0": ph0,
                        }
                        if key in params.dihedrals:
                            existing = params.dihedrals[key]
                            if isinstance(existing.get("multi"), list):
                                existing["multi"].append(entry)
                            else:
                                params.dihedrals[key] = {
                                    "multi": [existing, entry]
                                }
                        else:
                            params.dihedrals[key] = entry
                elif section == "IMPROPER":
                    cols = stripped.split()
                    if len(cols) >= 6 and not cols[0].startswith("!"):
                        t1, t2, t3, t4, kimp, psi0 = cols[0], cols[1], cols[2], cols[3], cols[4], cols[5]
                        key = (t1, t2, t3, t4)
                        kp = float(kimp)
                        ph0 = float(psi0)
                        if self.convert_units:
                            kp *= KJPERKCAL  # kcal/mol → kJ/mol
                            ph0 *= TORAD     # degrees → radians
                        params.impropers[key] = {
                            "k_phi": kp,
                            "phi0": ph0,
                        }
                elif section == "NONBONDED":
                    cols = stripped.split()
                    if len(cols) >= 3 and not cols[0].startswith("!"):
                        t1, eps, rmin = cols[0], cols[1], cols[2]
                        eps_val = abs(float(eps)) * KJPERKCAL  # kcal → kJ
                        sig_val = (float(rmin) * 2.0) / (10.0 * 1.122462048309373)
                        params.atom_types[t1] = {
                            "epsilon": eps_val,
                            "sigma": sig_val,
                        }
        return params


register("charmm", CharmmParser)
