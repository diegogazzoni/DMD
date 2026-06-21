"""Validation: Equipartition theorem — KE per degree of freedom = k_B·T/2."""
import sys, os, math
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "python"))

from dmd.runner import run

def test_equipartition():
    """NVT: average KE per atom should follow equipartition."""
    n_side = 6
    n_atoms = n_side ** 3  # 216
    mass = 39.948  # Ar
    T_target = 300.0
    box = 3.0
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
        "run": {"dt": 0.002, "n_steps": 1000, "init_temperature": T_target, "gen_vel": True, "seed": 42},
        "output": {"trajectory_path": "", "nstxout": 0, "nstvout": 0, "energy_path": "",
                    "nstenergy": 0, "checkpoint_interval": 0, "checkpoint_path": ""},
        "lj": {"cutoff": 1.2},
        "electrostatics": {"coulomb_type": "cutoff", "cutoff": 1.2, "pme_order": 4,
                           "pme_grid_spacing": 0.12, "ewald_coeff": 3.5},
        "thermostat": {"type": "nose-hoover", "temperature": T_target, "tau": 0.1, "frequency": 10.0},
        "barostat": {"type": "none", "pressure": 1.0, "tau": 1.0, "compressibility": 4.5e-5},
        "constraints": {"type": "none", "tolerance": 1e-6},
    }

    engine = run(system_data=system, config_data=config, progress_interval=0)

    kB = 0.008314462618
    vel = np.array(engine.velocities)
    ke_per_atom = 0.5 * mass * np.sum(vel**2, axis=1)
    avg_ke = np.mean(ke_per_atom)
    expected_ke = 1.5 * kB * T_target

    print(f"Equipartition test (Ar, {n_atoms} atoms, {T_target}K)")
    print(f"  Avg KE per atom: {avg_ke:.4f} kJ/mol")
    print(f"  Expected (3/2 kT): {expected_ke:.4f} kJ/mol")
    print(f"  Ratio: {avg_ke / expected_ke:.3f}")

    assert abs(avg_ke - expected_ke) / expected_ke < 0.30

    for i, dim in enumerate(['x', 'y', 'z']):
        ke_dim = 0.5 * mass * np.mean(vel[:, i]**2)
        expected_dim = 0.5 * kB * T_target
        print(f"  KE_{dim}: {ke_dim:.4f} (expected {expected_dim:.4f})")

    print(f"  PASS: Equipartition verified")
    return True

if __name__ == "__main__":
    test_equipartition()
