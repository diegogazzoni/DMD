import json
from _dmd_core import SimulationConfig

_SCHEMA = {
    "run": {
        "dt": float,
        "n_steps": int,
        "init_temperature": float,
        "gen_vel": bool,
        "seed": int,
    },
    "output": {
        "trajectory_path": str,
        "nstxout": int,
        "nstvout": int,
        "nstenergy": int,
        "energy_path": str,
        "checkpoint_interval": int,
        "checkpoint_path": str,
    },
    "lj": {
        "cutoff": float,
    },
    "electrostatics": {
        "coulomb_type": str,
        "cutoff": float,
        "pme_order": int,
        "pme_grid_spacing": float,
        "ewald_coeff": float,
    },
    "thermostat": {
        "type": str,
        "temperature": float,
        "tau": float,
        "frequency": float,
    },
    "barostat": {
        "type": str,
        "pressure": float,
        "tau": float,
        "compressibility": float,
    },
    "constraints": {
        "type": str,
        "tolerance": float,
    },
}

_JSON_TYPES = {
    float: "float",
    int: "integer",
    str: "string",
    bool: "boolean",
}


def _typename(t):
    return _JSON_TYPES.get(t, str(t))


def validate_config(config: dict) -> None:
    for sec_name, sec_schema in _SCHEMA.items():
        if sec_name not in config:
            raise ValueError(f"Missing required section '{sec_name}'")
        sec = config[sec_name]
        for key, expected_type in sec_schema.items():
            if key not in sec:
                raise ValueError(f"Missing required key '{sec_name}.{key}'")
            val = sec[key]
            if not isinstance(val, expected_type):
                raise ValueError(
                    f"Type mismatch for '{sec_name}.{key}': "
                    f"expected {_typename(expected_type)}, got {_typename(type(val))}"
                )
    for sec_name in config:
        if sec_name not in _SCHEMA:
            raise ValueError(f"Unknown section '{sec_name}'")
        for key in config[sec_name]:
            if key not in _SCHEMA[sec_name]:
                raise ValueError(f"Unknown key '{sec_name}.{key}'")


def apply_config(cfg: SimulationConfig, config: dict) -> None:
    validate_config(config)

    r = config["run"]
    cfg.dt = r["dt"]
    cfg.n_steps = r["n_steps"]
    cfg.init_temperature = r["init_temperature"]
    cfg.gen_vel = r["gen_vel"]
    cfg.seed = r["seed"]

    o = config["output"]
    cfg.trajectory_path = o["trajectory_path"]
    cfg.nstxout = o["nstxout"]
    cfg.nstvout = o["nstvout"]
    cfg.nstenergy = o["nstenergy"]
    cfg.energy_path = o["energy_path"]
    cfg.checkpoint_interval = o["checkpoint_interval"]
    cfg.checkpoint_path = o["checkpoint_path"]

    lj = config["lj"]
    cfg.lj_cutoff = lj["cutoff"]
    cfg.use_lj = True

    es = config["electrostatics"]
    cfg.coulomb_type = es["coulomb_type"]
    cfg.coulomb_cutoff = es["cutoff"]
    cfg.pme_order = es["pme_order"]
    cfg.pme_grid_spacing = es["pme_grid_spacing"]
    cfg.ewald_coeff = es["ewald_coeff"]

    tc = config["thermostat"]
    cfg.thermostat.type = tc["type"]
    cfg.thermostat.temperature = tc["temperature"]
    cfg.thermostat.tau = tc["tau"]
    cfg.thermostat.frequency = tc["frequency"]

    bc = config["barostat"]
    cfg.barostat.type = bc["type"]
    cfg.barostat.pressure = bc["pressure"]
    cfg.barostat.tau = bc["tau"]
    cfg.barostat.compressibility = bc["compressibility"]

    cc = config["constraints"]
    cfg.constraint_type = cc["type"]
    cfg.constraint_tolerance = cc["tolerance"]


def generate_config_template() -> str:
    t = {
        "run": {
            "dt": 0.002,
            "n_steps": 10000,
            "init_temperature": 300.0,
            "gen_vel": False,
            "seed": 42,
        },
        "output": {
            "trajectory_path": "trajectory.h5",
            "nstxout": 100,
            "nstvout": 100,
            "nstenergy": 10,
            "energy_path": "energy.h5",
            "checkpoint_interval": 500,
            "checkpoint_path": "checkpoint.bin",
        },
        "lj": {
            "cutoff": 1.2,
        },
        "electrostatics": {
            "coulomb_type": "cutoff",
            "cutoff": 1.2,
            "pme_order": 4,
            "pme_grid_spacing": 0.12,
            "ewald_coeff": 0.34,
        },
        "thermostat": {
            "type": "none",
            "temperature": 300.0,
            "tau": 0.1,
            "frequency": 10.0,
        },
        "barostat": {
            "type": "none",
            "pressure": 1.0,
            "tau": 1.0,
            "compressibility": 4.5e-5,
        },
        "constraints": {
            "type": "none",
            "tolerance": 1e-6,
        },
    }
    return json.dumps(t, indent=2)
