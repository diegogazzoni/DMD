import json
import numpy as np
from _dmd_core import Engine


def write_checkpoint(path: str, engine: Engine, system_data: dict) -> None:
    state = {
        "current_step": engine.current_step,
        "current_time": engine.current_time,
        "potential_energy": engine.potential_energy,
        "box_size": system_data["cell"]["box_size"],
        "positions": np.array(engine.positions).tolist(),
        "velocities": np.array(engine.velocities).tolist(),
    }
    payload = {
        "format": "dmd_checkpoint",
        "version": 1,
        "system": system_data,
        "state": state,
    }
    with open(path, "w") as f:
        json.dump(payload, f, indent=2)


def read_checkpoint(path: str) -> dict:
    with open(path) as f:
        data = json.load(f)
    if data.get("format") != "dmd_checkpoint":
        raise ValueError(f"Not a DMD checkpoint file: {path}")
    return data


def checkpoint_to_system_data(checkpoint: dict) -> dict:
    system = dict(checkpoint["system"])
    state = checkpoint["state"]
    system["cell"]["positions"] = state["positions"]
    system["cell"]["box_size"] = state["box_size"]
    return system
