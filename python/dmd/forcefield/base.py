from dataclasses import dataclass, field
from typing import Dict, Tuple


@dataclass
class FFParams:
    atom_types: Dict[str, dict] = field(default_factory=dict)
    bonds: Dict[Tuple[str, str], dict] = field(default_factory=dict)
    angles: Dict[Tuple[str, str, str], dict] = field(default_factory=dict)
    dihedrals: Dict[Tuple[str, str, str, str], dict] = field(default_factory=dict)
    impropers: Dict[Tuple[str, str, str, str], dict] = field(default_factory=dict)
