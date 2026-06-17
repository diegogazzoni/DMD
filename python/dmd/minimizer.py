import numpy as np
from _dmd_core import SimulationConfig, Engine
from dmd.system import SystemBuilder
from dmd.checkpoint import write_checkpoint


def minimize(
    system_data: dict,
    n_steps: int = 200,
    step_size: float = 0.001,
    energy_tol: float = 1e-4,
    output: str | None = None,
    verbose: bool = True,
) -> Engine:
    """Steepest descent energy minimization.

    Args:
        system_data: system.json dict
        n_steps: maximum number of minimization steps
        step_size: initial step size in nm
        energy_tol: stop if relative energy change < this
        output: optional checkpoint path to write final state
        verbose: print progress

    Returns:
        Engine after minimization (run with n_steps=1, forces at final positions)
    """
    pos = np.array(system_data["cell"]["positions"], dtype=np.float64)
    n_atoms = pos.shape[0]
    masses = np.array(system_data["atoms"]["mass"], dtype=np.float64)

    prev_pe = None

    for step in range(n_steps):
        cfg = _build_cfg(system_data, pos)
        engine = Engine(cfg)
        engine.run()

        pe = engine.potential_energy
        forces = np.array(engine.forces)

        if prev_pe is not None:
            delta = abs(pe - prev_pe) / max(abs(prev_pe), 1e-10)
            if delta < energy_tol:
                if verbose:
                    print(f"  Converged at step {step}: PE={pe:.6f}, Δ={delta:.2e}")
                break

        prev_pe = pe

        f_max = np.max(np.linalg.norm(forces, axis=1))
        if f_max < 1e-10:
            if verbose:
                print(f"  Zero forces at step {step}")
            break

        descent = forces / f_max
        pos = pos - step_size * descent

        if verbose and (step % 50 == 0 or step == n_steps - 1):
            print(f"  Step {step:4d}: PE={pe:.6f}, |F|_max={f_max:.4f}")

    cfg = _build_cfg(system_data, pos)
    engine = Engine(cfg)
    engine.run()

    if output:
        ckpt_system = dict(system_data)
        ckpt_system["cell"]["positions"] = pos.tolist()
        _write_minimizer_checkpoint(output, engine, ckpt_system)

    return engine


def _build_cfg(system_data: dict, positions: np.ndarray) -> SimulationConfig:
    data = dict(system_data)
    data["cell"]["positions"] = positions.tolist()
    builder = SystemBuilder(data)
    cfg = builder.build()
    cfg.gen_vel = False
    cfg.n_steps = 1
    return cfg


def _write_minimizer_checkpoint(path: str, engine: Engine, system_data: dict) -> None:
    import json
    import numpy as np
    state = {
        "current_step": 0,
        "current_time": 0.0,
        "potential_energy": engine.potential_energy,
        "box_size": system_data["cell"]["box_size"],
        "positions": np.array(system_data["cell"]["positions"]).tolist(),
        "velocities": [[0.0, 0.0, 0.0]] * len(system_data["cell"]["positions"]),
    }
    payload = {
        "format": "dmd_checkpoint",
        "version": 1,
        "system": system_data,
        "state": state,
    }
    with open(path, "w") as f:
        json.dump(payload, f, indent=2)
