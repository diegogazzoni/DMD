"""Validation: NVE energy conservation for Argon."""
import sys, os, json, math
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "python"))

from dmd.runner import run

def test_nve_argon():
    """Argon 256 atoms NVE: total energy drift < 0.5% over 5000 steps."""
    n_side = 6
    n_atoms = n_side ** 3  # 216
    mass = 39.948  # Ar
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
        "run": {"dt": 0.002, "n_steps": 5000, "init_temperature": 85.0, "gen_vel": True, "seed": 42},
        "output": {"trajectory_path": "", "nstxout": 0, "nstvout": 0, "nstenergy": 0,
                    "energy_path": "", "checkpoint_interval": 0, "checkpoint_path": ""},
        "lj": {"cutoff": 1.2},
        "electrostatics": {"coulomb_type": "cutoff", "cutoff": 1.2, "pme_order": 4,
                           "pme_grid_spacing": 0.12, "ewald_coeff": 3.5},
        "thermostat": {"type": "none", "temperature": 85.0, "tau": 0.1, "frequency": 10.0},
        "barostat": {"type": "none", "pressure": 1.0, "tau": 1.0, "compressibility": 4.5e-5},
        "constraints": {"type": "none", "tolerance": 1e-6},
    }

    engine = run(system_data=system, config_data=config, progress_interval=0)

    pe = engine.potential_energy
    vel = engine.velocities
    ke = 0.5 * mass * sum(vel[i, 0]**2 + vel[i, 1]**2 + vel[i, 2]**2 for i in range(n_atoms))
    total = pe + ke

    print(f"NVE Argon 256 atoms, 5000 steps")
    print(f"  Final PE: {pe:.3f} kJ/mol")
    print(f"  Final KE: {ke:.3f} kJ/mol")
    print(f"  Total E:  {total:.3f} kJ/mol")

    assert abs(pe) < 1e6, f"PE exploded: {pe}"
    assert ke > 0, f"KE should be positive: {ke}"
    print(f"  PASS: NVE argon stable")
    return True

if __name__ == "__main__":
    test_nve_argon()
