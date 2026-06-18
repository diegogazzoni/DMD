# DMD — Modular Molecular Dynamics Engine

DMD is a CPU-based molecular dynamics engine written in modern C++20 with a
Python API. It runs NVE, NVT, and NPT simulations with support for multiple
force fields, thermostats, and barostats. Input is provided as structured
JSON (system.json + config.json); trajectory output uses the HDF5-based
H5MD standard.

## Features

- **Integrator:** Velocity Verlet (symplectic, time-reversible)
- **Force fields:**
  - Lennard-Jones (with multi-type support and Verlet neighbor list)
  - Coulomb (direct Ewald summation and Smooth Particle Mesh Ewald)
  - Harmonic bonds, angles, and periodic dihedrals
- **Thermostats:** Berendsen (velocity rescaling), Nose-Hoover (extended Lagrangian), Andersen (stochastic collisions)
- **Barostats:** Berendsen (isotropic scaling), Andersen (extended Lagrangian piston)
- **Input:** Structured JSON (system.json + config.json)
- **Output:** H5MD trajectory files (HDF5-based, standard format for MD analysis)
- **Checkpoint/restart:** Bit-exact restart for reproducible trajectories
- **Configuration:** Strict JSON schema (no implicit defaults, all keys required)
- **Python API:** SystemBuilder, format converters (PDB, PSF, GRO), force field parsers (CHARMM)

## Architecture

```
src/
├── sim/            SimulationEngine, SimulationConfig, checkpoint
├── force/          force components: LJ, CoulombDirect, CoulombPME, bonds, angles, dihedrals
├── integrate/      Velocity Verlet integrator
├── thermostat/     Berendsen, Nose-Hoover, Andersen
├── barostat/       Berendsen, Andersen
├── trajectory/     H5MD writer
└── core/           Cell, SystemData, constants

python/
└── dmd/            Python API (SystemBuilder, converters, force field parsers, runner)
```

## Dependencies

| Dependency | Version | Notes |
|---|---|---|
| C++ compiler | C++20 | Clang 14+, GCC 12+, MSVC 2022 |
| CMake | ≥ 3.25 | Build system |
| FFTW3 | ≥ 3.3 | Required for Coulomb PME |
| HDF5 | ≥ 1.12 | Required for H5MD trajectory output |
| Google Test | 1.15.2 | Fetched automatically for tests |
| Python | ≥ 3.10 | Required for the Python API |

### macOS

```bash
# Install system dependencies via Homebrew
brew install cmake fftw hdf5

# Build
cmake -B build
cmake --build build -j$(sysctl -n hw.logicalcpu)
```

### Linux (Ubuntu / Debian)

```bash
# Install system dependencies
sudo apt install cmake g++ libfftw3-dev libhdf5-dev

# Build
cmake -B build
cmake --build build -j$(nproc)
```

### Linux (Fedora / RHEL)

```bash
sudo dnf install cmake gcc-c++ fftw-devel hdf5-devel
cmake -B build
cmake --build build -j$(nproc)
```

If HDF5 is installed in a non-standard location, set `CMAKE_PREFIX_PATH`:

```bash
cmake -B build -DCMAKE_PREFIX_PATH=/path/to/hdf5
```

## Quick Start

```python
import dmd, json

with open("system.json") as f:
    system = json.load(f)
with open("config.json") as f:
    config = json.load(f)

engine = dmd.run(system_data=system, config_data=config)
print(f"Potential energy: {engine.potential_energy:.4f} kJ/mol")
```

See the [User Guide](docs/user_guide.md) for the full tutorial.

## File Formats

### system.json

Human-readable JSON describing the atomic system: box size, positions, masses,
charges, atom types, force field parameters (LJ pairs, bonds, angles, dihedrals),
and exclusion pairs.

### config.json

Strict JSON schema with seven required sections. Every key must be present —
no implicit defaults. Use `dmd.generate_config_template()` to generate a
complete template.

### H5MD trajectory

Output follows the H5MD standard (HDF5-based):

```
/h5md/version                    [1, 1]
/particles/atoms/position/value  [n_frames, n_atoms, 3]
/particles/atoms/position/time   [n_frames]
/particles/atoms/position/step   [n_frames]
/particles/atoms/velocity/value  [n_frames, n_atoms, 3]  (if nstvout > 0)
/particles/atoms/box/edges/value [n_frames, 3, 3]
```

## Testing

```bash
ctest --test-dir build
```

DMD has **15+ unit and integration tests** covering all force components,
thermostats, barostats, integrator, checkpoint/restart, H5MD I/O.

## License

MIT
