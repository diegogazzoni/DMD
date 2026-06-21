# ADR 0003: Replace binary .dmdin with human-readable JSON

**Status:** Accepted (implemented 2026-06-18)

**Context:** Originally DMD used a custom binary format (`.dmdin`) for system input and a C++ JSON parser (nlohmann) for runtime config. This created a three-layer problem: system.dmdin → config.json → C++ internal state.

**Decision:** Eliminate all binary formats and C++ JSON parsing. The Python layer is the sole human interface:
- `system.json` replaces `.dmdin` (human-readable, diffable, version-controllable)
- `config.json` remains for runtime parameters but validation moves to pure Python (`config.py`)
- C++ core is JSON-agnostic — it receives a `SimulationConfig` struct via pybind11

**Consequences:**
- Removed: `src/sysbin/`, `src/main/`, `src/sim/json_config.*`, nlohmann_json dependency
- Added: `python/dmd/config.py` with `validate_config()`, `apply_config()`, `generate_config_template()`
- All input validation is now in Python (faster iteration, richer error messages)
- C++ build time reduced (~2500 lines removed, no nlohmann_json FetchContent)
