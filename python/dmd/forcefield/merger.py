"""Merge PSF connectivity with force field parameters."""

from typing import Dict, List, Optional, Tuple
from dmd.forcefield.base import FFParams


def _wildcard_match(key: Tuple[str, ...], param_key: Tuple[str, ...]) -> bool:
    """Match a key like ('CA','CB') against a param key like ('CA','X'),
    where 'X' is a CHARMM wildcard matching any type."""
    if len(key) != len(param_key):
        return False
    for a, b in zip(key, param_key):
        if a != b and b != "X":
            return False
    return True


def _best_match(
    key: Tuple[str, ...],
    param_dict: Dict,
) -> Optional[dict]:
    """Find the best matching param entry for a given type tuple.
    Prefers exact match, then wildcard."""
    if key in param_dict:
        entry = param_dict[key]
        if isinstance(entry, dict) and "multi" in entry:
            return entry["multi"]
        return entry
    # Try reverse if 2-tuple
    if len(key) == 2:
        rev = (key[1], key[0])
        if rev in param_dict:
            entry = param_dict[rev]
            if isinstance(entry, dict) and "multi" in entry:
                return entry["multi"]
            return entry
    # Try wildcard
    for pk, entry in param_dict.items():
        if _wildcard_match(key, pk):
            if isinstance(entry, dict) and "multi" in entry:
                return entry["multi"]
            return entry
        if len(key) == 2:
            if _wildcard_match((key[1], key[0]), pk):
                if isinstance(entry, dict) and "multi" in entry:
                    return entry["multi"]
                return entry
    return None


def merge_ff(
    psf_data: dict,
    ff_params: FFParams,
    atom_types: Optional[List[str]] = None,
) -> dict:
    """Merge PSF connectivity with FF parameters to produce a complete
    force_field section ready for system.json.

    Args:
        psf_data: output of from_psf() (contains force_field, optionally atom_types)
        ff_params: output of load_forcefield()
        atom_types: optional override — list of type names per atom index.
                    If not provided, taken from psf_data['atom_types'].

    Returns:
        dict with keys: bonds, angles, dihedrals, impropers (all with params)
    """
    if atom_types is None:
        atom_types = psf_data.get("atom_types", [])
    force_field = {}
    psf_ff = psf_data.get("force_field", {})

    # Bonds
    bonds_out = []
    for b in psf_ff.get("bonds", []):
        entry = dict(b)
        if atom_types:
            ti = atom_types[b["i"]]
            tj = atom_types[b["j"]]
            p = _best_match((ti, tj), ff_params.bonds)
            if p:
                entry.update(p)
        bonds_out.append(entry)
    if bonds_out:
        force_field["bonds"] = bonds_out

    # Angles
    angles_out = []
    for a in psf_ff.get("angles", []):
        entry = dict(a)
        if atom_types:
            ti = atom_types[a["i"]]
            tj = atom_types[a["j"]]
            tk = atom_types[a["k"]]
            p = _best_match((ti, tj, tk), ff_params.angles)
            if p:
                entry.update(p)
        angles_out.append(entry)
    if angles_out:
        force_field["angles"] = angles_out

    # Dihedrals
    dih_out = []
    for d in psf_ff.get("dihedrals", []):
        base = {"i": d["i"], "j": d["j"], "k": d["k"], "l": d["l"]}
        if atom_types:
            ti = atom_types[d["i"]]
            tj = atom_types[d["j"]]
            tk = atom_types[d["k"]]
            tl = atom_types[d["l"]]
            p = _best_match((ti, tj, tk, tl), ff_params.dihedrals)
            if isinstance(p, list):
                for term in p:
                    entry = dict(base)
                    entry.update(term)
                    dih_out.append(entry)
            elif p:
                entry = dict(base)
                entry.update(p)
                dih_out.append(entry)
        else:
            dih_out.append(base)
    if dih_out:
        force_field["dihedrals"] = dih_out

    # Impropers
    imp_out = []
    for im in psf_ff.get("impropers", []):
        entry = dict(im)
        if atom_types:
            ti = atom_types[im["i"]]
            tj = atom_types[im["j"]]
            tk = atom_types[im["k"]]
            tl = atom_types[im["l"]]
            p = _best_match((ti, tj, tk, tl), ff_params.impropers)
            if p:
                entry.update(p)
        imp_out.append(entry)
    if imp_out:
        force_field["impropers"] = imp_out

    # LJ from atom_types
    if atom_types and ff_params.atom_types:
        unique_types = sorted(set(atom_types))
        type_set = list(unique_types)
        n = len(type_set)
        pairs = []
        for i, ti in enumerate(type_set):
            for j, tj in enumerate(type_set):
                p_i = ff_params.atom_types.get(ti, {})
                p_j = ff_params.atom_types.get(tj, {})
                eps_i = p_i.get("epsilon", 0.0)
                eps_j = p_j.get("epsilon", 0.0)
                sig_i = p_i.get("sigma", 1.0)
                sig_j = p_j.get("sigma", 1.0)
                eps = (eps_i * eps_j) ** 0.5 if eps_i != 0.0 and eps_j != 0.0 else 0.0
                sig = (sig_i + sig_j) * 0.5
                pairs.append({
                    "type_i": ti,
                    "type_j": tj,
                    "sigma": sig,
                    "epsilon": eps,
                })
        force_field["lennard_jones"] = {"pairs": pairs}

    return {"force_field": force_field}


def generate_constraints(
    psf_data: dict,
    ff_params: FFParams,
    atom_types: Optional[List[str]] = None,
) -> List[dict]:
    """Generate SHAKE constraint list from PSF connectivity + FF parameters.

    For each water-like angle (i-j-k) with a central atom j, produces three
    constraints: i-j (bond), j-k (bond), i-k (angle-derived).

    Args:
        psf_data: output of from_psf()
        ff_params: output of load_forcefield()
        atom_types: optional override for atom type names

    Returns:
        list of dicts with keys: i, j, distance
    """
    if atom_types is None:
        atom_types = psf_data.get("atom_types", [])
    psf_ff = psf_data.get("force_field", {})

    # Build bond distance lookup: (sorted atom indices) -> distance
    bond_dist = {}
    for b in psf_ff.get("bonds", []):
        ai, aj = b["i"], b["j"]
        d = None
        if atom_types:
            ti = atom_types[ai]
            tj = atom_types[aj]
            p = _best_match((ti, tj), ff_params.bonds)
            if p and "r0" in p:
                d = p["r0"]
            elif p and "distance" in p:
                d = p["distance"]
        if d is None:
            d = b.get("r0", b.get("distance", 0.0))
        bond_dist[(min(ai, aj), max(ai, aj))] = d

    constraints = []
    seen_pairs = set()

    def add_constraint(i, j, distance):
        key = (min(i, j), max(i, j))
        if key not in seen_pairs:
            constraints.append({"i": i, "j": j, "distance": distance})
            seen_pairs.add(key)

    for angle in psf_ff.get("angles", []):
        ai, aj, ak = angle["i"], angle["j"], angle["k"]

        # Bond i-j
        d_ij = bond_dist.get((min(ai, aj), max(ai, aj)))
        if d_ij is not None:
            add_constraint(ai, aj, d_ij)

        # Bond j-k
        d_jk = bond_dist.get((min(aj, ak), max(aj, ak)))
        if d_jk is not None:
            add_constraint(aj, ak, d_jk)

        # Angle-derived i-k distance
        theta0 = None
        if atom_types:
            ti = atom_types[ai]
            tj = atom_types[aj]
            tk = atom_types[ak]
            p = _best_match((ti, tj, tk), ff_params.angles)
            if p and "theta0" in p:
                theta0 = p["theta0"]

        if theta0 is not None and d_ij is not None and d_jk is not None:
            d_ik = (d_ij * d_ij + d_jk * d_jk
                    - 2.0 * d_ij * d_jk * _cos(theta0)) ** 0.5
            add_constraint(ai, ak, d_ik)

    return constraints


def _cos(angle_rad: float) -> float:
    """Cosine, safe for small angles."""
    import math
    return math.cos(angle_rad)
