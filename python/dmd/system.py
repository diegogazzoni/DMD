from __future__ import annotations
import json
import numpy as np
from _dmd_core import SimulationConfig
from dmd.config import apply_config


class SystemBuilder:
    """Build a SimulationConfig from system.json and config.json inputs."""

    def __init__(self, system_json: dict | str | None = None):
        self.cfg = SimulationConfig()
        if system_json is not None:
            if isinstance(system_json, str):
                system_json = json.loads(system_json)
            self.from_system_json(system_json)

    def from_system_json(self, data: dict) -> "SystemBuilder":
        n = len(data["cell"]["positions"])
        self.cfg.box_size = data["cell"].get("box_size", 3.0)

        pos = np.array(data["cell"]["positions"], dtype=np.float64)
        self.cfg.pos_x = pos[:, 0].copy()
        self.cfg.pos_y = pos[:, 1].copy()
        self.cfg.pos_z = pos[:, 2].copy()

        masses = np.array(data["atoms"]["mass"], dtype=np.float64)
        self.cfg.mass = masses

        if "charge" in data["atoms"]:
            charges = np.array(data["atoms"]["charge"], dtype=np.float64)
        else:
            charges = np.zeros(n, dtype=np.float64)
        self.cfg.charges = charges

        type_names = data["atoms"].get("type", [])
        name_to_idx = {name: i for i, name in enumerate(sorted(set(type_names)))}
        atom_types = np.array([name_to_idx[t] for t in type_names], dtype=np.int32)
        self.cfg.atom_types = atom_types

        ff = data.get("force_field", {})
        lj = ff.get("lennard_jones", {})
        if lj:
            self.cfg.use_lj = True
            self.cfg.lj_cutoff = lj.get("cutoff", 2.5)
            n_types = len(set(type_names))
            sigma = np.full((n_types, n_types), 1.0, dtype=np.float64)
            epsilon = np.full((n_types, n_types), 0.0, dtype=np.float64)
            for p in lj.get("pairs", []):
                i = name_to_idx[p["type_i"]]
                j = name_to_idx[p["type_j"]]
                sigma[i][j] = sigma[j][i] = p.get("sigma", 1.0)
                epsilon[i][j] = epsilon[j][i] = p.get("epsilon", 0.0)
            self.cfg.lj_sigma = sigma.flatten()
            self.cfg.lj_epsilon = epsilon.flatten()

        bonds = ff.get("bonds", [])
        if bonds:
            self.cfg.bond_i = np.array([b["i"] for b in bonds], dtype=np.int32)
            self.cfg.bond_j = np.array([b["j"] for b in bonds], dtype=np.int32)
            if "k" in bonds[0]:
                self.cfg.bond_k = np.array([b.get("k", 0.0) for b in bonds], dtype=np.float64)
            if "r0" in bonds[0]:
                self.cfg.bond_r0 = np.array([b.get("r0", 0.0) for b in bonds], dtype=np.float64)

        angles = ff.get("angles", [])
        if angles:
            self.cfg.angle_i = np.array([a["i"] for a in angles], dtype=np.int32)
            self.cfg.angle_j = np.array([a["j"] for a in angles], dtype=np.int32)
            self.cfg.angle_k = np.array([a["k"] for a in angles], dtype=np.int32)
            if "k_theta" in angles[0]:
                self.cfg.angle_k_theta = np.array([a.get("k_theta", 0.0) for a in angles], dtype=np.float64)
            if "theta0" in angles[0]:
                self.cfg.angle_theta0 = np.array([a.get("theta0", 0.0) for a in angles], dtype=np.float64)

        dihedrals = ff.get("dihedrals", [])
        if dihedrals:
            self.cfg.dihedral_i = np.array([d["i"] for d in dihedrals], dtype=np.int32)
            self.cfg.dihedral_j = np.array([d["j"] for d in dihedrals], dtype=np.int32)
            self.cfg.dihedral_k = np.array([d["k"] for d in dihedrals], dtype=np.int32)
            self.cfg.dihedral_l = np.array([d["l"] for d in dihedrals], dtype=np.int32)
            if "k_phi" in dihedrals[0]:
                self.cfg.dihedral_k_phi = np.array([d.get("k_phi", 0.0) for d in dihedrals], dtype=np.float64)
            if "periodicity" in dihedrals[0]:
                self.cfg.dihedral_periodicity = np.array([d.get("periodicity", 1) for d in dihedrals], dtype=np.int32)
            if "phi0" in dihedrals[0]:
                self.cfg.dihedral_phi0 = np.array([d.get("phi0", 0.0) for d in dihedrals], dtype=np.float64)

        impropers = ff.get("impropers", [])
        if impropers:
            self.cfg.improper_i = np.array([im["i"] for im in impropers], dtype=np.int32)
            self.cfg.improper_j = np.array([im["j"] for im in impropers], dtype=np.int32)
            self.cfg.improper_k = np.array([im["k"] for im in impropers], dtype=np.int32)
            self.cfg.improper_l = np.array([im["l"] for im in impropers], dtype=np.int32)
            if "k_phi" in impropers[0]:
                self.cfg.improper_k_phi = np.array([im.get("k_phi", 0.0) for im in impropers], dtype=np.float64)
            if "phi0" in impropers[0]:
                self.cfg.improper_phi0 = np.array([im.get("phi0", 0.0) for im in impropers], dtype=np.float64)

        constraints = ff.get("constraints", [])
        if constraints:
            self.cfg.constraint_i = [c["i"] for c in constraints]
            self.cfg.constraint_j = [c["j"] for c in constraints]
            self.cfg.constraint_dist = [c["distance"] for c in constraints]

        exclusions = data.get("exclusions", [])
        if exclusions:
            self.cfg.excl_i = [e[0] for e in exclusions]
            self.cfg.excl_j = [e[1] for e in exclusions]
        else:
            self.cfg.excl_i = []
            self.cfg.excl_j = []

        return self

    def apply_config_json(self, path_or_dict: str | dict) -> "SystemBuilder":
        if isinstance(path_or_dict, dict):
            apply_config(self.cfg, path_or_dict)
        else:
            with open(path_or_dict) as f:
                apply_config(self.cfg, json.load(f))
        return self

    def build(self) -> SimulationConfig:
        return self.cfg
