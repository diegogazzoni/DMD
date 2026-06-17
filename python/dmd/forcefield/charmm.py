"""CHARMM PAR/PRM force field file parser."""

import re
from typing import Dict
from dmd.forcefield.base import FFParams
from dmd.forcefield.registry import register


class CharmmParser:
    """Parse CHARMM PAR/PRM force field parameter files."""

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
                        params.bonds[(t1, t2)] = params.bonds[(t2, t1)] = {
                            "k": float(kb),
                            "r0": float(b0),
                        }
                elif section == "ANGLES":
                    cols = stripped.split()
                    if len(cols) >= 6 and not cols[0].startswith("!"):
                        t1, t2, t3, ktheta, theta0 = cols[0], cols[1], cols[2], cols[3], cols[4]
                        key = (t1, t2, t3)
                        params.angles[key] = {
                            "k_theta": float(ktheta),
                            "theta0": float(theta0),
                        }
                elif section == "DIHEDRALS":
                    cols = stripped.split()
                    if len(cols) >= 7 and not cols[0].startswith("!"):
                        t1, t2, t3, t4 = cols[0], cols[1], cols[2], cols[3]
                        kchi, periodicity, delta = cols[4], cols[5], cols[6]
                        key = (t1, t2, t3, t4)
                        entry = {
                            "k_phi": float(kchi),
                            "periodicity": int(periodicity),
                            "phi0": float(delta),
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
                        params.impropers[key] = {
                            "k_phi": float(kimp),
                            "phi0": float(psi0),
                        }
                elif section == "NONBONDED":
                    cols = stripped.split()
                    if len(cols) >= 3 and not cols[0].startswith("!"):
                        t1, eps, rmin = cols[0], cols[1], cols[2]
                        params.atom_types[t1] = {
                            "epsilon": float(eps),
                            "sigma": float(rmin) * 2.0,
                        }
        return params


register("charmm", CharmmParser)
