import json
import os
import tempfile
from _dmd_core import (
    SimulationConfig,
    Engine,
    read_dmdin,
    write_dmdin,
    apply_json_config,
)
from dmd.system import SystemBuilder


def run(
    system_data: str | dict | None = None,
    config_data: str | dict | None = None,
    system_json: str | None = None,
    config_json: str | None = None,
    dmdin_path: str | None = None,
    output_dmdin: str | None = None,
    output_traj: str | None = None,
) -> Engine:
    """Orchestrate a DMD simulation run.

    Provide one of:
        - system_json + config_json (file paths)
        - system_data + config_data (dicts or JSON strings)
        - dmdin_path (prebuilt .dmdin file + optional config_json)
    """
    if dmdin_path:
        cfg = read_dmdin(dmdin_path)
        if config_json:
            apply_json_config(cfg, config_json)
        elif config_data:
            with tempfile.NamedTemporaryFile(
                mode="w", suffix=".json", delete=False
            ) as f:
                if isinstance(config_data, dict):
                    json.dump(config_data, f)
                else:
                    f.write(config_data)
                tmpcfg = f.name
            try:
                apply_json_config(cfg, tmpcfg)
            finally:
                os.unlink(tmpcfg)
    else:
        if isinstance(system_data, str):
            system_data = json.loads(system_data)
        builder = SystemBuilder(system_data)
        if config_data:
            if isinstance(config_data, str):
                config_data = json.loads(config_data)
            builder.apply_config_json(config_data)
        cfg = builder.build()

    if output_dmdin:
        write_dmdin(output_dmdin, cfg)

    engine = Engine(cfg)
    engine.run()

    if output_traj:
        _write_traj(engine, output_traj)

    return engine


def _write_traj(engine: Engine, path: str) -> None:
    import numpy as np
    pos = np.array(engine.positions)
    np.save(path, pos)
