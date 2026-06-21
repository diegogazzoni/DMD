"""Validation: TIP3P water geometry and property tests."""
import sys, os, json, math
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "python"))

from dmd.convert import from_pdb, from_psf
from dmd.forcefield import load_forcefield, merge_ff, generate_constraints
from dmd.runner import run

BASE = os.path.join(os.path.dirname(__file__), "..", "..", "test_acqua")

def test_water_geometry():
    """SHAKE must maintain TIP3P geometry throughout simulation."""
    pdb = from_pdb(os.path.join(BASE, "water_box.pdb"))
    psf = from_psf(os.path.join(BASE, "water_box.psf"), scale_14=1.0)
    ff = load_forcefield(os.path.join(BASE, "tip3p.prm"), format="charmm")

    positions = np.array(pdb["cell"]["positions"], dtype=np.float64) / 10.0
    system = {
        "cell": {"positions": positions.tolist(), "box_size": pdb["cell"]["box_size"] / 10.0},
        "atoms": {"mass": pdb["atoms"]["mass"], "charge": pdb["atoms"]["charge"], "type": psf["atom_types"]},
    }
    constraints = generate_constraints(psf, ff)
    system["force_field"] = {"constraints": constraints}
    merged = merge_ff(psf, ff)
    system["force_field"]["lennard_jones"] = merged["force_field"]["lennard_jones"]
    if "exclusions" in psf:
        system["exclusions"] = psf["exclusions"]

    config_data = json.load(open(os.path.join(BASE, "config.json")))

    engine = run(system_data=system, config_data=config_data, progress_interval=0)

    n_atoms = len(pdb["cell"]["positions"])
    pos = np.array(engine.positions)
    n_waters = n_atoms // 3

    # Measure all waters at final frame
    oh1_dists = np.linalg.norm(pos[1::3] - pos[0::3], axis=1)
    oh2_dists = np.linalg.norm(pos[2::3] - pos[0::3], axis=1)
    hh_dists = np.linalg.norm(pos[2::3] - pos[1::3], axis=1)

    cos_theta = (np.sum((pos[1::3] - pos[0::3]) * (pos[2::3] - pos[0::3]), axis=1)
                 / (oh1_dists * oh2_dists))
    angles = np.degrees(np.arccos(cos_theta))

    mean_oh = np.mean(oh1_dists)
    std_oh = np.std(oh1_dists)
    mean_hh = np.mean(hh_dists)
    mean_angle = np.mean(angles)

    print(f"TIP3P water geometry ({n_waters} molecules, {config_data['run']['n_steps']} steps)")
    print(f"  O-H distance: {mean_oh:.5f} ± {std_oh:.5f} nm (target 0.09572)")
    print(f"  H-H distance: {mean_hh:.5f} nm (target 0.15139)")
    print(f"  H-O-H angle:  {mean_angle:.2f}° (target 104.52°)")
    print(f"  PE: {engine.potential_energy:.1f} kJ/mol")

    # Assertions
    assert abs(mean_oh - 0.09572) < 1e-4, f"O-H distance {mean_oh:.5f} != 0.09572"
    assert abs(mean_hh - 0.15139) < 1e-4, f"H-H distance {mean_hh:.5f} != 0.15139"
    assert abs(mean_angle - 104.52) < 0.5, f"Angle {mean_angle:.2f} != 104.52"

    print(f"  PASS: Water geometry verified through {config_data['run']['n_steps']} steps")
    return True

def test_pme_energy_stability():
    """PME water should have stable potential energy (~1 kJ/mol/atom)."""
    pdb = from_pdb(os.path.join(BASE, "water_box.pdb"))
    psf = from_psf(os.path.join(BASE, "water_box.psf"), scale_14=1.0)
    ff = load_forcefield(os.path.join(BASE, "tip3p.prm"), format="charmm")

    positions = np.array(pdb["cell"]["positions"], dtype=np.float64) / 10.0
    system = {
        "cell": {"positions": positions.tolist(), "box_size": pdb["cell"]["box_size"] / 10.0},
        "atoms": {"mass": pdb["atoms"]["mass"], "charge": pdb["atoms"]["charge"], "type": psf["atom_types"]},
    }
    constraints = generate_constraints(psf, ff)
    system["force_field"] = {"constraints": constraints}
    merged = merge_ff(psf, ff)
    system["force_field"]["lennard_jones"] = merged["force_field"]["lennard_jones"]
    if "exclusions" in psf:
        system["exclusions"] = psf["exclusions"]

    config_data = json.load(open(os.path.join(BASE, "config.json")))

    engine = run(system_data=system, config_data=config_data, progress_interval=0)

    n_atoms = len(pdb["cell"]["positions"])
    pe = engine.potential_energy
    pe_per_atom = pe / n_atoms
    forces = np.array(engine.forces)
    max_force = np.max(np.abs(forces))

    print(f"PME water energy stability ({n_atoms} atoms)")
    print(f"  PE: {pe:.1f} kJ/mol ({pe_per_atom:.2f} kJ/mol/atom)")
    print(f"  Max |force|: {max_force:.1f} kJ/(mol·nm)")

    # Reasonable bounds for TIP3P water at 300K
    assert 0 < pe_per_atom < 5.0, f"PE per atom {pe_per_atom:.2f} out of range [0, 5] kJ/mol"
    assert 0 < max_force < 10000, f"Max force {max_force:.0f} too high"

    print(f"  PASS: Energy stable")
    return True

if __name__ == "__main__":
    test_water_geometry()
    test_pme_energy_stability()
