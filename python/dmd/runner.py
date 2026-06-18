from __future__ import annotations
import json
from _dmd_core import SimulationConfig, Engine
from dmd.system import SystemBuilder
from dmd.checkpoint import read_checkpoint, checkpoint_to_system_data
from dmd.status import SimulationStatus
from dmd.config import apply_config


def run(
    system_data: str | dict | None = None,
    config_data: str | dict | None = None,
    system_json: str | None = None,
    config_json: str | None = None,
    checkpoint: str | None = None,
    progress_interval: int = 100,
) -> Engine:
    """Orchestrate a DMD simulation run.

    Provide one of:
        - system_json + config_json (file paths)
        - system_data + config_data (dicts or JSON strings)
        - checkpoint (continuation from a checkpoint.json file)

    Args:
        progress_interval: print status every N steps (0 to disable)
    """
    if checkpoint:
        ckpt = read_checkpoint(checkpoint)
        system_data = checkpoint_to_system_data(ckpt)

    if system_data:
        if isinstance(system_data, str):
            system_data = json.loads(system_data)
        builder = SystemBuilder(system_data)
        if config_data:
            if isinstance(config_data, str):
                config_data = json.loads(config_data)
            builder.apply_config_json(config_data)
        cfg = builder.build()
    else:
        raise ValueError("Provide one of: system_data or checkpoint")

    engine = Engine(cfg)

    if progress_interval > 0 and cfg.n_steps > 0:
        status = SimulationStatus(engine, dt=cfg.dt, interval=progress_interval)
        status.start()
        chunk = progress_interval
        remaining = cfg.n_steps
        while remaining > 0:
            take = min(chunk, remaining)
            engine.run_n(take)
            remaining -= take
            status.step()
        status.finish()
    else:
        engine.run()

    return engine
