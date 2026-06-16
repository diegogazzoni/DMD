# DMD — Modular Molecular Dynamics Engine

DMD is a CPU-based molecular dynamics engine written in modern C++20. It runs NVE, NVT, and NPT simulations with support for multiple force fields, thermostats, and barostats. Input is provided as a compact binary format (`.dmdin`) paired with a strict JSON configuration file; trajectory output uses the HDF5-based H5MD standard.

## Features

- **Integrator:** Velocity Verlet (symplectic, time-reversible)
- **Force fields:**
  - Lennard-Jones (with multi-type support and Verlet neighbor list)
  - Coulomb (direct Ewald summation and Smooth Particle Mesh Ewald)
  - Harmonic bonds, angles, and periodic dihedrals
- **Thermostats:** Berendsen (velocity rescaling), Nose-Hoover (extended Lagrangian), Andersen (stochastic collisions)
- **Barostats:** Berendsen (isotropic scaling), Andersen (extended Lagrangian piston)
- **Input:** Binary `.dmdin` format (compact, self-describing, no parsing ambiguity)
- **Output:** H5MD trajectory files (HDF5-based, standard format for MD analysis)
- **Checkpoint/restart:** Bit-exact restart for reproducible trajectories
- **Configuration:** Strict JSON schema (no implicit defaults, all keys required)
- **CLI:** `dmd run` and `dmd config --template`

## Architecture

DMD is organised in five layers, each depending only on the layer below:

```
CLI / Tests
      |
Simulation Engine      run loop, checkpoint, trajectory I/O
      |
Force Engine           force components (LJ, Coulomb, bonded)
      |
Thermostat / Barostat  coupling algorithms
      |
System Layer           Cell, SystemData, binary I/O
```

```
src/
├── main/           CLI entry point
├── sim/            SimulationEngine, SimulationConfig, checkpoint, JSON config
├── force/          force components: LJ, CoulombDirect, CoulombPME, bonds, angles, dihedrals
├── integrate/      Velocity Verlet integrator
├── thermostat/     Berendsen, Nose-Hoover, Andersen
├── barostat/       Berendsen, Andersen
├── trajectory/     H5MD writer
├── sysbin/         .dmdin binary reader/writer
└── core/           Cell, SystemData, constants
```

## Dependencies

| Dependency | Version | Notes |
|---|---|---|
| C++ compiler | C++20 | Clang 14+, GCC 12+, MSVC 2022 |
| CMake | ≥ 3.25 | Build system |
| FFTW3 | ≥ 3.3 | Required for Coulomb PME |
| HDF5 | ≥ 1.12 | Required for H5MD trajectory output |
| nlohmann-json | 3.11.3 | Fetched automatically by CMake |
| Google Test | 1.15.2 | Fetched automatically for tests |

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
# Install system dependencies
sudo dnf install cmake gcc-c++ fftw-devel hdf5-devel

# Build
cmake -B build
cmake --build build -j$(nproc)
```

If HDF5 is installed in a non-standard location, set `CMAKE_PREFIX_PATH`:

```bash
cmake -B build -DCMAKE_PREFIX_PATH=/path/to/hdf5
```

## Quick Start

```bash
# 1. Generate a config.json template
./build/src/main/dmd config --template > config.json

# 2. Edit config.json with your parameters
#    (all keys are required — replace example values)

# 3. Run a simulation
./build/src/main/dmd run system.dmdin config.json
```

To run the tutorial (Argon NVE), see the [User Guide](docs/user_guide.md).

## File Formats

### .dmdin (binary system input)

A compact binary format containing everything needed to describe a system:

- **Header** (120 bytes): magic (`0x444D444E`), version, atom count, cell type, number of atom types, position flags, FF section offset, position section offset, box matrix (9 doubles)
- **System section**: per-type masses and charges
- **FF section**: bonded parameters (bonds, angles, dihedrals) and LJ cross-type parameters
- **Positions section**: positions and velocities (stored as 32-bit floats for compactness) and per-atom types (32-bit ints)

The binary format eliminates parsing ambiguity, is faster to read than text formats, and produces smaller files.

### config.json (run configuration)

Strict JSON schema with seven required sections. Every key must be present — no implicit defaults. Use `dmd config --template` to generate a complete template.

Full documentation in the [User Guide](docs/user_guide.md#configuration-reference).

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

## CLI Reference

```
dmd run <system.dmdin> <config.json>
    Run a simulation from a system file and configuration.

dmd config --template
    Print a complete config.json template to stdout.

dmd --help
    Show usage information.
```

## Testing

```bash
ctest --test-dir build
```

DMD has **20+ unit and integration tests** covering all force components, thermostats, barostats, integrator, checkpoint/restart, H5MD I/O, and the full production pipeline.

## License

MIT
