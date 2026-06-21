from typing import Dict, Type
from dmd.forcefield.base import FFParams


_registry: Dict[str, type] = {}


def register(name: str, parser_cls: type):
    _registry[name] = parser_cls


def load_forcefield(path: str, format: str = "charmm", **kwargs) -> FFParams:
    if format not in _registry:
        available = ", ".join(_registry.keys())
        raise ValueError(
            f"Unknown force field format '{format}'. "
            f"Available: {available or 'none'}"
        )
    parser = _registry[format](**kwargs)
    return parser.parse(path)
