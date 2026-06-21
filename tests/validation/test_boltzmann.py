"""Validation: Maxwell-Boltzmann velocity distribution in NVT."""
import sys, os, json, math
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "python"))

from dmd.runner import run

def test_boltzmann_distribution():
    """NVT: velocity distribution should follow Maxwell-Boltzmann."""
    n_side = 8
    n_atoms = n_side ** 3  # 512
    mass = 39.948  # Ar
    T_target = 300.0
    box = 4.0
    L = box / 2.0
    spacing = box / n_side
    pos = [[i * spacing - L, j * spacing - L, k * spacing - L]
           for i in range(n_side) for j in range(n_side) for k in range(n_side)]

    system = {
        "cell": {"positions": pos, "box_size": box},
        "atoms": {"mass": [mass] * n_atoms, "charge": [0.0] * n_atoms, "type": ["Ar"] * n_atoms},
        "force_field": {
            "lennard_jones": {
                "pairs": [{"type_i": "Ar", "type_j": "Ar", "sigma": 0.3405, "epsilon": 0.996}],
                "cutoff": 1.2,
            }
        },
    }

    config = {
        "run": {"dt": 0.002, "n_steps": 2000, "init_temperature": T_target, "gen_vel": True, "seed": 7},
        "output": {"trajectory_path": "", "nstxout": 0, "nstvout": 200, "energy_path": "",
                    "nstenergy": 0, "checkpoint_interval": 0, "checkpoint_path": ""},
        "lj": {"cutoff": 1.2},
        "electrostatics": {"coulomb_type": "cutoff", "cutoff": 1.2, "pme_order": 4,
                           "pme_grid_spacing": 0.12, "ewald_coeff": 3.5},
        "thermostat": {"type": "nose-hoover", "temperature": T_target, "tau": 0.1, "frequency": 10.0},
        "barostat": {"type": "none", "pressure": 1.0, "tau": 1.0, "compressibility": 4.5e-5},
        "constraints": {"type": "none", "tolerance": 1e-6},
    }

    engine = run(system_data=system, config_data=config, progress_interval=0)

    vel = np.array(engine.velocities)
    n = vel.shape[0]
    kB = 0.008314462618

    # Compute temperature
    ke_total = 0.5 * mass * np.sum(vel**2)
    T_inst = 2.0 * ke_total / (3.0 * n) / kB

    sigma_v = np.sqrt(kB * T_target / mass)
    v_all = vel.flatten()
    v_std = np.std(v_all)

    print(f"Boltzmann test (Ar, {n} atoms, {T_target}K, 2000 steps)")
    print(f"  Measured T: {T_inst:.1f} K (target {T_target} K)")
    print(f"  Velocity std: {v_std:.4f} (expected {sigma_v:.4f})")
    print(f"  Velocity mean: {np.mean(v_all):.6f} (expected ~0)")

    assert abs(T_inst - T_target) / T_target < 0.3
    assert abs(v_std - sigma_v) / sigma_v < 0.2
    assert abs(np.mean(v_all)) < 0.01

    print(f"  PASS: Boltzmann distribution verified")
    return True

if __name__ == "__main__":
    test_boltzmann_distribution()
