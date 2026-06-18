# DMD — Complete User Guide

## 1. Introduction to Molecular Dynamics

### 1.1 What is Molecular Dynamics?

**Molecular Dynamics** (MD) is a computer simulation technique that studies the
motion of atoms and molecules over time. Each atom is treated as a point particle
with a mass and charge; its motion obeys classical mechanics (Newton's second law):

```
F = m · a
```

By numerically integrating this equation for every atom at each time step
(Δt ~ 1–2 femtoseconds), we obtain a **trajectory**: the evolution of the system
over nanoseconds or microseconds.

MD can compute thermodynamic (temperature, pressure, energy), structural (radial
distribution functions, bond angles) and dynamical (diffusion coefficients,
vibrational spectra) properties of systems ranging from a few hundred to millions
of atoms.

### 1.2 Potentials and Force Fields

Forces between atoms are derived from **force fields**, semi-empirical potential
functions that approximate the system's potential energy as a sum of contributions:

```
E_tot = E_bond + E_angle + E_dihedral + E_LJ + E_Coulomb
```

| Contribution | Description | Model |
|---|---|---|
| **E_bond** | Chemical bond between two atoms | Harmonic: ½k(r − r₀)² |
| **E_angle** | Angle between three atoms | Harmonic: ½k_θ(θ − θ₀)² |
| **E_dihedral** | Rotation around a bond | Periodic: k_φ[1+cos(nφ−φ₀)] |
| **E_LJ** (Lennard-Jones) | van der Waals forces | 4ε[(σ/r)¹² − (σ/r)⁶] |
| **E_Coulomb** | Electrostatic interactions | k_e · q_i · q_j / r |

**Bond/angle/dihedral** are called **bonded forces** because they involve atoms
connected by covalent bonds. **LJ and Coulomb** are called **non-bonded forces**
because they act between all atom pairs regardless of connectivity.

### 1.3 Short-range vs Long-range

The most important computational distinction is between **short-range** and
**long-range** forces:

- **Short-range** (LJ, direct Coulomb): decay rapidly with r. A **cutoff**
  is applied (typically 1.0–1.4 nm): beyond that distance the contribution is
  negligible. Complexity is O(N · M) where M is the average number of neighbors
  within the cutoff (~50–500 atoms).

- **Long-range** (Coulomb): the electrostatic interaction decays as 1/r, too
  slowly to be truncated. Using a cutoff causes severe artifacts. Methods like
  **Ewald summation** or **Particle Mesh Ewald (PME)** split the potential into
  a real-space part (computed with a cutoff) and a reciprocal part (computed in
  Fourier space via 3D FFT). Complexity: O(N log N).

### 1.4 Constraints

In classical MD, bonds involving hydrogens (C–H, O–H, N–H) vibrate at very high
frequencies (~3000 cm⁻¹). Integrating them correctly would require Δt < 0.5 fs,
slowing the simulation. The solution is to **constrain** these bond lengths to
fixed values using algorithms like **SHAKE** or **LINCS**, allowing Δt = 2 fs
even in the presence of hydrogens.

Constraints are applied _after_ integration by correcting positions and/or velocities.

### 1.5 Thermodynamic Ensembles

| Ensemble | Constants | Description |
|---|---|---|
| **NVE** | Number of atoms, Volume, Energy | Microcanonical: isolated system, total energy conserved |
| **NVT** | Number of atoms, Volume, Temperature | Canonical: constant temperature via thermostat |
| **NPT** | Number of atoms, Pressure, Temperature | Isobaric-isothermal: constant pressure and temperature |

NVT requires a **thermostat** that exchanges heat with a virtual heat bath.
NPT adds a **barostat** that scales the volume to maintain the target pressure.

---

## 2. Software Architecture

DMD is organized in two layers:

```
┌─────────────────────────────────────────┐
│            Python API Layer              │  system.py, runner.py, minimizer.py
│  convert/ (PDB, PSF, GRO)               │  checkpoint.py, forcefield/ (CHARMM)
│  forcefield/ (parser, merger)           │
├─────────────────────────────────────────┤
│   pybind11 bridge (core.cpp)            │  _dmd_core.so
├─────────────────────────────────────────┤
│   SimulationEngine                      │  sim/simulation_engine.cpp
│     ForceEngine                         │  force/force_engine.cpp
│       HarmonicBond                      │  force/harmonic_bond.cpp
│       HarmonicAngle                     │  force/harmonic_angle.cpp
│       PeriodicDihedral                  │  force/periodic_dihedral.cpp
│       LennardJones                      │  force/lennard_jones.cpp
│       CoulombDirect                     │  force/coulomb_direct.cpp
│       CoulombPME                        │  force/coulomb_pme.cpp
│       CoulombExclusion                  │  force/coulomb_exclusion.cpp
│     Integrator (Velocity Verlet)        │  integrate/integrator.cpp
│     Thermostat (Berendsen/NH/Andersen)  │  thermostat/*.cpp
│     Barostat (Berendsen/Andersen)       │  barostat/*.cpp
│     H5MDWriter                          │  trajectory/h5md_writer.cpp
│     Checkpoint (binary)                 │  sim/checkpoint.cpp
├─────────────────────────────────────────┤
│   Cell / SystemData / Constants         │  core/*.cpp
├─────────────────────────────────────────┤
│            C++ Core Layer               │
└─────────────────────────────────────────┘
```

### Architectural Principles

- **Clean separation**: C++ handles the hot loop (forces, integration, trajectory
  I/O). Python handles setup, format conversion, orchestration.
- **SystemData SoA**: all atomic properties are in parallel vectors (Structure of
  Arrays): `pos_x`, `pos_y`, `pos_z` instead of `pos[i].x`. Maximizes cache locality
  and facilitates span-based access.
- **ForceComponent**: each energy term is a class implementing
  `compute(span pos_x/y/z, cell, span forces_x/y/z, energy&)`. Registered at
  startup in a fixed order.
- **Strict JSON config**: the config.json file follows a rigid schema: every key
  must be present, unknown keys produce errors. No hidden defaults.
- **Native H5MD**: trajectories are in H5MD (HDF5), an open standard. No
  proprietary formats.

### Units and Conversions

All internal quantities use the **DMD** system (coherent with GROMACS):
- Length: **nm**
- Energy: **kJ/mol**
- Force: **kJ/(mol·nm)**
- Time: **ps**
- Charge: **elementary charge e⁻**

Force field parsers (e.g. CHARMM, which uses kcal/mol and Å) convert
automatically on import.

---

## 3. The system.json Format

The **system.json** file describes the atomic system: the simulation box, atomic
positions, masses, charges, atom types, and force field parameters. It can be
created manually or generated by converters (PDB, PSF, GRO).

### 3.1 General Structure

```json
{
    "cell": {
        "positions": [[x1, y1, z1], [x2, y2, z2], ...],
        "box_size": 3.0
    },
    "atoms": {
        "mass": [39.948, 39.948, ...],
        "charge": [0.0, 0.0, ...],
        "type": ["Ar", "Ar", ...]
    },
    "force_field": {
        "lennard_jones": { ... },
        "bonds": [ ... ],
        "angles": [ ... ],
        "dihedrals": [ ... ],
        "impropers": [ ... ]
    },
    "exclusions": [[0, 1], [0, 2], [1, 3], ...]
}
```

### 3.2 `cell` Section

| Field | Type | Required | Description |
|---|---|---|---|
| `positions` | array of [x,y,z] | yes | Atomic coordinates in nm |
| `box_size` | float | yes | Side length of the cubic box in nm |

The box is currently **cubic** (single side length). Triclinic, orthorhombic,
and rhombic dodecahedron cells are planned for future versions.

### 3.3 `atoms` Section

| Field | Type | Required | Description |
|---|---|---|---|
| `mass` | float array | yes | Atomic masses in g/mol |
| `charge` | float array | yes | Atomic charges in e⁻ |
| `type` | string array | yes | Atom type names (e.g. "OH2", "HW1", "Ar") |

All three arrays must have the same length (N atoms). Type names must match
those used in `force_field.lennard_jones.pairs`.

### 3.4 `force_field` Section

#### Lennard-Jones

```json
"lennard_jones": {
    "pairs": [
        {"type_i": "Ar", "type_j": "Ar", "sigma": 0.3405, "epsilon": 0.997},
        {"type_i": "OH2", "type_j": "OH2", "sigma": 0.3150, "epsilon": 0.636}
    ]
}
```

| Field | Type | Required | Description |
|---|---|---|---|
| `type_i`, `type_j` | string | yes | Atom type names |
| `sigma` | float | yes | LJ σ parameter in **nm** |
| `epsilon` | float | yes | LJ ε parameter in **kJ/mol** |

Important: σ is the true LJ parameter (potential minimum at 2^(1/6)·σ),
**not** Rmin/2 (the van der Waals radius). When using CHARMM parameters,
the `charmm.py` parser converts automatically:
- Rmin/2 (Å) → σ (nm): `σ = (2 · Rmin/2) / (10 · 2^(1/6))`
- ε (kcal/mol) → ε (kJ/mol): `ε × 4.184`

Combination rules for mixed pairs (if not explicitly specified) follow
Lorentz-Berthelot:
- `σ_ij = (σ_i + σ_j) / 2`  (arithmetic mean)
- `ε_ij = √(ε_i · ε_j)`     (geometric mean)

Pairs are expanded into an `n_types × n_types` matrix at startup.

#### Bonds

```json
"bonds": [
    {"i": 0, "j": 1, "k": 502416.0, "r0": 0.09572}
]
```

| Field | Type | Required | Description |
|---|---|---|---|
| `i`, `j` | int | yes | Atom indices (0-based) |
| `k` | float | yes | Force constant in kJ/(mol·nm²) |
| `r0` | float | yes | Equilibrium length in nm |

Energy: `E = ½ · k · (r − r₀)²` where `r = |r_i − r_j|`.

#### Angles

```json
"angles": [
    {"i": 0, "j": 1, "k": 2, "k_theta": 628.02, "theta0": 1.824}
]
```

| Field | Type | Required | Description |
|---|---|---|---|
| `i`, `j`, `k` | int | yes | Atom indices, j = central atom |
| `k_theta` | float | yes | Force constant in kJ/(mol·rad²) |
| `theta0` | float | yes | Equilibrium angle in **radians** |

Energy: `E = ½ · k_θ · (θ − θ₀)²`.

#### Proper Dihedrals

```json
"dihedrals": [
    {"i": 0, "j": 1, "k": 2, "l": 3, "k_phi": 5.0, "periodicity": 3, "phi0": 0.0}
]
```

| Field | Type | Required | Description |
|---|---|---|---|
| `i`, `j`, `k`, `l` | int | yes | Atom indices |
| `k_phi` | float | yes | Force constant in kJ/mol |
| `periodicity` | int | yes | Multiplicity n |
| `phi0` | float | yes | Phase in **radians** |

Energy: `E = k_φ · [1 + cos(n · φ − φ₀)]`.

Multi-term dihedrals: multiple entries with the same quadruple (i,j,k,l) but
different `k_phi`, `periodicity`, `phi0`. The CHARMM parser handles these
automatically.

#### Impropers

```json
"impropers": [
    {"i": 0, "j": 1, "k": 2, "l": 3, "k_phi": 10.0, "phi0": 3.14159}
]
```

Same fields as proper dihedrals. They describe out-of-plane torsions
(e.g. a planar central atom in a carbonyl group). Computed by the same
`PeriodicDihedral` component.

### 3.5 `exclusions` Section

```json
"exclusions": [[0, 1], [0, 2], [1, 3]]
```

Atom pairs for which **non-bonded** interactions (LJ + Coulomb) are **not**
computed. Typically:
- **1-2 exclusion**: directly bonded atoms
- **1-3 exclusion**: atoms separated by two bonds (angle)
- **1-4 exclusion** (optional): scales dihedral interactions

The PSF parser (`from_psf()`) automatically generates exclusions from
connectivity (bond → 1-2, angle → 1-3) and can include 1-4 scaling via
`scale_14`. PDB and GRO parsers do not generate exclusions.

### 3.6 Generating system.json

#### From PDB + PSF + Force Field (complete workflow)

```python
import dmd
import json

# 1. Read coordinates
system = dmd.from_pdb("structure.pdb")

# 2. Read topology (bonds, angles, dihedrals, exclusions)
psf = dmd.from_psf("structure.psf", scale_14=1.0)
system["atoms"]["type"] = psf["atom_types"]
system["force_field"].update(psf["force_field"])

# 3. Load force field parameters
ff_params = dmd.load_forcefield("params.prm", format="charmm")

# 4. Merge: match parameters to PSF atom types
merged = dmd.merge_ff(system, ff_params)
system["force_field"].update(merged["force_field"])
system["exclusions"] = psf["exclusions"]

# 5. Save
with open("system.json", "w") as f:
    json.dump(system, f, indent=2)
```

#### From GRO (coordinates only, no bonded terms)

```python
system = dmd.from_gro("structure.gro")
```

#### From scratch (manual, for testing)

```python
system = {
    "cell": {
        "positions": [[0.0, 0.0, 0.0], [0.375, 0.375, 0.0], ...],
        "box_size": 1.5
    },
    "atoms": {
        "mass": [39.948] * 32,
        "charge": [0.0] * 32,
        "type": ["Ar"] * 32
    },
    "force_field": {
        "lennard_jones": {
            "pairs": [
                {"type_i": "Ar", "type_j": "Ar", "sigma": 0.3405, "epsilon": 0.997}
            ]
        }
    }
}
```

---

---

## 4. The config.json Format — All Parameters

The **config.json** file contains ALL simulation parameters. The schema is
**strict**: every key is mandatory. To generate a template:
```python
import dmd
print(dmd.generate_config_template())
```

### 5.1 `run` Section — Execution Parameters

```json
"run": {
    "dt": 0.001,
    "n_steps": 5000,
    "init_temperature": 300.0,
    "gen_vel": true,
    "seed": 42
}
```

| Key | Type | Required | Description |
|---|---|---|---|
| `dt` | float | yes | Integration time step in **picoseconds**. Typical values: 0.002 ps (2 fs) for systems with constrained C–H bonds; 0.001 ps for all-atom systems without constraints. |
| `n_steps` | int | yes | Total number of integration steps. Total simulated time = dt × n_steps. |
| `init_temperature` | float | yes | Velocity initialization temperature in **Kelvin**. Used only if `gen_vel = true`. |
| `gen_vel` | bool | yes | If `true`, generates initial velocities from a Maxwell-Boltzmann distribution at `init_temperature`. If `false`, velocities start at zero (useful for minimization or restart). |
| `seed` | int | yes | Random number generator seed. Same seed = same initial velocity sequence. Changing the seed yields different initial conditions. |

### 5.2 `output` Section — Output and Trajectory

```json
"output": {
    "trajectory_path": "trajectory.h5",
    "nstxout": 100,
    "nstvout": 0,
    "nstenergy": 10,
    "energy_path": "energy.h5",
    "checkpoint_interval": 500,
    "checkpoint_path": "checkpoint.bin"
}
```

| Key | Type | Required | Description |
|---|---|---|---|
| `trajectory_path` | string | yes | Path to H5MD trajectory file. Empty string `""` = no trajectory output. |
| `nstxout` | int | yes | Step interval for writing **positions** to the H5MD file. 0 = don't write positions. Typical: 50–1000. |
| `nstvout` | int | yes | Step interval for writing **velocities** to the H5MD file. 0 = don't write velocities (saves disk space). |
| `nstenergy` | int | yes | Step interval for writing **energies** to a separate H5MD file. 0 = don't write. |
| `energy_path` | string | yes | Path to H5MD energy file. Used if `nstenergy > 0`. |
| `checkpoint_interval` | int | yes | Step interval between checkpoints (state saves). 0 = disabled. |
| `checkpoint_path` | string | yes | Path to the binary checkpoint file. Used if `checkpoint_interval > 0` and for the final checkpoint. |

**Relationship between nstxout and nstvout**: positions and velocities are
written at **the same frames** (same step). If `nstxout = 50` and
`nstvout = 50`, each frame contains both. If `nstvout = 0`, velocities
are not written. Different frequencies are not supported — keep them equal.

### 5.3 `lj` Section — Lennard-Jones

```json
"lj": {
    "cutoff": 1.2
}
```

| Key | Type | Required | Description |
|---|---|---|---|
| `cutoff` | float | yes | Lennard-Jones cutoff radius in **nm**. Typical: 1.0–1.4 nm. Beyond this distance the LJ contribution is considered zero. |

The Verlet list (neighbor list within cutoff + skin) is rebuilt every 10
integration steps.

### 5.4 `electrostatics` Section

```json
"electrostatics": {
    "coulomb_type": "pme",
    "cutoff": 1.2,
    "pme_order": 4,
    "pme_grid_spacing": 0.12,
    "ewald_coeff": 3.5
}
```

| Key | Type | Required | Description |
|---|---|---|---|
| `coulomb_type` | string | yes | Electrostatic treatment. Values: `"cutoff"` (simple truncation), `"direct"` (Ewald direct sum with erfc), `"pme"` (full Particle Mesh Ewald). |
| `cutoff` | float | yes | Real-space Ewald cutoff radius in nm. Must be ≥ the LJ cutoff. |
| `pme_order` | int | yes | B-spline order for PME grid interpolation. Typical: 4 (cubic) or 5 (more accurate, more expensive). |
| `pme_grid_spacing` | float | yes | Desired PME grid spacing in nm. Actual grid points are computed as `ceil(L / spacing)` where L is the box size. Typical: 0.10–0.15 nm. |
| `ewald_coeff` | float | yes | Ewald split coefficient α in **nm⁻¹**. Typical values: 3.0–4.0. Large α → more real-space work, less reciprocal work. Small α → more FFT work. The optimal choice balances the two costs. |

**Choosing `coulomb_type`:**
- `"cutoff"`: hard truncation at `cutoff`. **Warning**: produces artifacts for
  charged systems. Use only for testing or neutral systems without charges.
- `"direct"`: Ewald direct sum with erfc/r term. Includes the exclusion
  correction (1-2, 1-3). More accurate than cutoff but O(N²) complexity.
- `"pme"`: Particle Mesh Ewald. Splits the potential into real-space (erfc/r
  within cutoff) + reciprocal (3D FFT on a grid) + self-energy.
  Complexity O(N log N). **Recommended** for systems > 2000 atoms.

**The Ewald coefficient α** determines how much of the force is computed in
real vs reciprocal space. A reasonable value is `α = 3.5 / cutoff`.

### 5.5 `thermostat` Section

```json
"thermostat": {
    "type": "nose-hoover",
    "temperature": 300.0,
    "tau": 0.1,
    "frequency": 10.0
}
```

| Key | Type | Required | Description |
|---|---|---|---|
| `type` | string | yes | Thermostat type. Values: `"none"`, `"berendsen"`, `"nose-hoover"`, `"andersen"`. |
| `temperature` | float | yes | Target temperature in **Kelvin**. |
| `tau` | float | yes | Coupling constant in **ps** for Berendsen and Nose-Hoover. Smaller = stronger coupling. |
| `frequency` | float | yes | Collision frequency in **ps⁻¹** for Andersen. Average number of collisions per atom per picosecond. |

**Thermostat details:**

| Type | Mechanism | Exact ensemble? | Typical use |
|---|---|---|---|
| `"none"` | No thermostat (NVE) | Yes (microcanonical) | Energy conservation tests, isolated systems |
| `"berendsen"` | Proportional velocity rescaling | No | Rapid equilibration (small τ = 0.1 ps), **never in production** |
| `"nose-hoover"` | Extended variable ζ (friction) | Yes (canonical NVT) | Production (τ = 0.5–2.0 ps) |
| `"andersen"` | Stochastic collisions with heat bath | Yes (canonical NVT) | Small systems, when discontinuous dynamics is acceptable |

**Nose-Hoover**: the exact canonical thermostat recommended for production.
It introduces an additional degree of freedom ζ (friction coefficient) with
mass `Q = N_f · k_B · T₀ · τ²`. The dynamics is deterministic and reversible.

**Berendsen**: scales velocities as `λ = √(1 + Δt/τ · (T₀/T − 1))`.
Does not produce a true canonical ensemble (suppresses temperature fluctuations).
Use only for initial equilibration.

**Andersen**: at each step, every atom has probability `P = ν · Δt` of
undergoing a collision that replaces its velocity with a Maxwell-Boltzmann
sample. Produces a true canonical ensemble but disrupts dynamics (velocities
are not continuous).

### 5.6 `barostat` Section

```json
"barostat": {
    "type": "none",
    "pressure": 1.0,
    "tau": 1.0,
    "compressibility": 4.5e-5
}
```

| Key | Type | Required | Description |
|---|---|---|---|
| `type` | string | yes | Barostat type. Values: `"none"`, `"berendsen"`, `"andersen"`. |
| `pressure` | float | yes | Target pressure in **bar**. |
| `tau` | float | yes | Barostat coupling constant in **ps**. |
| `compressibility` | float | yes | Isothermal compressibility in **bar⁻¹** (Berendsen only). For water at 300K: 4.5×10⁻⁵ bar⁻¹. |

**Barostat details:**

| Type | Mechanism | Typical use |
|---|---|---|
| `"none"` | Fixed volume (NVE or NVT) | Constant volume simulations |
| `"berendsen"` | Proportional volume scaling | Equilibration, **not exact** |
| `"andersen"` | Extended piston (ϵ degree of freedom) | NPT production |

**Berendsen barostat**: scales volume as `μ = 1 − Δt/τ_P · β · (P₀ − P)`.
Position and box scale factor: `μ^(1/3)`. Does not reproduce correct volume
fluctuations in the NPT ensemble.

**Andersen barostat**: uses an extended Lagrangian with piston variable ϵ
of mass `W = N_f · k_B · T₀ · τ²`. Positions scale as `exp(ϵ̇ · Δt)`.
More accurate than Berendsen for the NPT ensemble.

### 5.7 `constraints` Section

```json
"constraints": {
    "type": "none",
    "tolerance": 1e-6
}
```

| Key | Type | Required | Description |
|---|---|---|---|
| `type` | string | yes | Constraint type. Currently only `"none"`. |
| `tolerance` | float | yes | Constraint solver tolerance (used by SHAKE). |

SHAKE constraints are implemented in C++ but not yet connected to the
SimulationEngine. Hydrogen bonds must currently be integrated with a small
time step (dt = 0.001 ps).

### 5.8 Complete Examples

#### NVE (microcanonical)

```json
{
    "run": {"dt": 0.002, "n_steps": 5000, "init_temperature": 85.0, "gen_vel": true, "seed": 42},
    "output": {"trajectory_path": "", "nstxout": 0, "nstvout": 0, "nstenergy": 0, "energy_path": "",
               "checkpoint_interval": 0, "checkpoint_path": ""},
    "lj": {"cutoff": 1.2},
    "electrostatics": {"coulomb_type": "cutoff", "cutoff": 1.2, "pme_order": 4, "pme_grid_spacing": 0.12, "ewald_coeff": 3.5},
    "thermostat": {"type": "none", "temperature": 85.0, "tau": 0.1, "frequency": 10.0},
    "barostat": {"type": "none", "pressure": 1.0, "tau": 1.0, "compressibility": 4.5e-5},
    "constraints": {"type": "none", "tolerance": 1e-6}
}
```

#### NVT with Nose-Hoover (production)

```json
{
    "run": {"dt": 0.001, "n_steps": 50000, "init_temperature": 300.0, "gen_vel": true, "seed": 42},
    "output": {"trajectory_path": "traj.h5", "nstxout": 100, "nstvout": 100, "nstenergy": 10,
               "energy_path": "energy.h5", "checkpoint_interval": 1000, "checkpoint_path": "checkpoint.bin"},
    "lj": {"cutoff": 1.2},
    "electrostatics": {"coulomb_type": "pme", "cutoff": 1.2, "pme_order": 4, "pme_grid_spacing": 0.12, "ewald_coeff": 3.5},
    "thermostat": {"type": "nose-hoover", "temperature": 300.0, "tau": 1.0, "frequency": 10.0},
    "barostat": {"type": "none", "pressure": 1.0, "tau": 1.0, "compressibility": 4.5e-5},
    "constraints": {"type": "none", "tolerance": 1e-6}
}
```

#### NPT with Berendsen (equilibration)

```json
{
    "run": {"dt": 0.002, "n_steps": 10000, "init_temperature": 120.0, "gen_vel": true, "seed": 42},
    "output": {"trajectory_path": "", "nstxout": 0, "nstvout": 0, "nstenergy": 10, "energy_path": "",
               "checkpoint_interval": 500, "checkpoint_path": "checkpoint.bin"},
    "lj": {"cutoff": 1.2},
    "electrostatics": {"coulomb_type": "cutoff", "cutoff": 1.2, "pme_order": 4, "pme_grid_spacing": 0.12, "ewald_coeff": 3.5},
    "thermostat": {"type": "berendsen", "temperature": 120.0, "tau": 0.1, "frequency": 10.0},
    "barostat": {"type": "berendsen", "pressure": 100.0, "tau": 1.0, "compressibility": 4.5e-5},
    "constraints": {"type": "none", "tolerance": 1e-6}
}
```

---

## 5. Running a Simulation — End-to-End

### 6.1 Complete Workflow

```
┌─────────────────┐
│   system.json   │  ← from PDB+PSF+FF, or GRO, or hand-written
└────────┬────────┘
         ↓
┌─────────────────┐
│  dmd.run()      │  ← reads system.json + config.json
│  SystemBuilder  │     builds SimulationConfig
└────────┬────────┘
         ↓
┌─────────────────┐
│  Engine.run()   │  ← MD loop: forces → integration → thermostat → output
│                  │     writes H5MD + checkpoint
└────────┬────────┘
         ↓
┌─────────────────┐
│  Results        │  → engine.positions, .velocities, .potential_energy
│  H5MD file      │  → trajectory for analysis (MDAnalysis, VMD, Python)
│  Checkpoint     │  → restart
└─────────────────┘
```

### 6.2 Minimal Example Script

```python
import dmd

with open("system.json") as f:
    import json
    system = json.load(f)

with open("config.json") as f:
    config = json.load(f)

engine = dmd.run(system_data=system, config_data=config)

print(f"Steps: {engine.current_step}")
print(f"Time: {engine.current_time:.4f} ps")
print(f"Potential energy: {engine.potential_energy:.4f} kJ/mol")
```

### 6.3 Output with Progress Bar

```python
engine = dmd.run(
    system_data=system,
    config_data=config,
    progress_interval=100
)
```

Output format:
```
        Step  Time (ps)    PE (kJ/mol)   Box (nm)   ns/day
--------------------------------------------------------
        500     1.0000       -85.2341     3.0000 175575.5
       1000     2.0000       -85.1234     3.0000 180234.1
```

### 6.4 Complete PDB+PSF+FF Script

```python
import dmd, json

# 1. Prepare system
system = dmd.from_pdb("structure.pdb")

psf = dmd.from_psf("structure.psf", scale_14=1.0)
system["atoms"]["type"] = psf["atom_types"]
system["force_field"].update(psf["force_field"])

ff_params = dmd.load_forcefield("params.prm", format="charmm")
merged = dmd.merge_ff(system, ff_params)
system["force_field"].update(merged["force_field"])
system["exclusions"] = psf["exclusions"]

# 2. Configure simulation
config = {
    "run": {"dt": 0.001, "n_steps": 10000, "init_temperature": 300.0,
            "gen_vel": True, "seed": 42},
    "output": {"trajectory_path": "traj.h5", "nstxout": 50, "nstvout": 0,
               "nstenergy": 10, "energy_path": "energy.h5",
               "checkpoint_interval": 1000, "checkpoint_path": "checkpoint.bin"},
    "lj": {"cutoff": 1.2},
    "electrostatics": {"coulomb_type": "pme", "cutoff": 1.2,
                       "pme_order": 4, "pme_grid_spacing": 0.12, "ewald_coeff": 3.5},
    "thermostat": {"type": "nose-hoover", "temperature": 300.0,
                   "tau": 1.0, "frequency": 10.0},
    "barostat": {"type": "none", "pressure": 1.0, "tau": 1.0,
                 "compressibility": 4.5e-5},
    "constraints": {"type": "none", "tolerance": 1e-6},
}

# 3. Run
engine = dmd.run(system_data=system, config_data=config, progress_interval=100)

# 4. Basic analysis
import numpy as np
pos = np.array(engine.positions)
vel = np.array(engine.velocities)
masses = np.array(system["atoms"]["mass"])
ekin = 0.5 * np.sum(masses * np.sum(vel**2, axis=1))
kb = 0.008314462618
T = 2.0 * ekin / (3.0 * len(masses) * kb)
print(f"Final temperature: {T:.2f} K")
```

### 5.5 Custom analysis

```python
import numpy as np

# After engine.run():
pos = np.array(engine.positions)
vel = np.array(engine.velocities)
masses = np.array(system["atoms"]["mass"])

# Kinetic energy
ekin = 0.5 * np.sum(masses * np.sum(vel**2, axis=1))

# Instantaneous temperature
kb = 0.008314462618  # kJ/(mol·K)
T = 2.0 * ekin / (3.0 * len(masses) * kb)
print(f"Temperature: {T:.2f} K")
```

### 5.6 Energy Minimization

Before a production simulation, it is good practice to **minimize the energy**
to remove unfavorable steric contacts:

```python
import dmd

# Minimize
info = dmd.minimize(
    system_data=system,
    n_steps=500,
    step_size=0.001,
    energy_tol=0.01,
    output={"checkpoint_path": "minimized_checkpoint.json"}
)

# Start simulation from the minimum
engine = dmd.run(
    checkpoint="minimized_checkpoint.json",
    config_data=config
)
```

The algorithm is **steepest descent with backtracking**: if energy increases,
the step is reduced (×0.5); if it decreases, the step is increased (×1.2).
Converges when the energy change falls below `energy_tol` or `n_steps` is reached.

### 5.7 Example: TIP3P Water

The complete example for a TIP3P water box (522 molecules, 1566 atoms) is in
`test_acqua/`:

```bash
cd test_acqua

# Generate system (PDB + PSF + PRM)
python3 generate_water_box.py

# Run simulation
python3 run_water_test.py

# Visualize trajectory
python3 visualize_vdw.py
```

TIP3P parameters used:
- **O–O LJ**: σ = 0.3150 nm, ε = 0.636 kJ/mol
- **O–H bond**: k = 502416 kJ/(mol·nm²), r₀ = 0.09572 nm
- **H–O–H angle**: k_θ = 628.02 kJ/(mol·rad²), θ₀ = 104.52°
- **Charges**: O = −0.834 e⁻, H = +0.417 e⁻
- **Box**: 2.5 nm cubic, 522 H₂O

---

## 6. Algorithms in Detail

### 6.1 Velocity Verlet Integrator

DMD uses the **Velocity Verlet** integrator, a symplectic (energy-conserving)
and time-reversible algorithm:

```
v(t + ½Δt) = v(t) + (Δt/2) · a(t)              # half-kick
r(t + Δt)  = r(t) + Δt · v(t + ½Δt)            # advance
a(t + Δt)  = F(t + Δt) / m                     # forces
v(t + Δt)  = v(t + ½Δt) + (Δt/2) · a(t + Δt)  # half-kick
```

Properties:
- **Symplectic**: conserves total energy in the limit of infinitesimal steps.
  At finite steps, energy oscillates around a constant value without drift.
- **Reversible**: flipping the sign of Δt moves backward in time.
- **Local error**: O(Δt³) on positions, O(Δt²) on velocities.
- **Global error**: O(Δt²) on positions.

The full MD loop in DMD:

```
for step in range(n_steps):
    integrator.half_kick(dt/2)              # v ← v + dt/2 · F/m
    integrator.advance(dt, box)             # r ← r + dt · v (with PBC)
    force_engine.compute(pos → F, E)        # compute all forces
    integrator.half_kick(dt/2)              # v ← v + dt/2 · F/m
    thermostat.apply(vel, dt)               # regulate temperature
    barostat.apply(pos, box, dt)            # regulate pressure
    constraints.apply(pos, vel)             # constraints (SHAKE)
    trajectory.write(pos, vel, box, step)   # H5MD
    checkpoint.save(state)                  # periodic checkpoint
```

### 6.2 Nose-Hoover Thermostat

The **Nose-Hoover** thermostat extends the system with an additional degree of
freedom ζ (friction coefficient). The equations of motion are:

```
ṙ_i = v_i
v̇_i = F_i/m_i − ζ · v_i
ζ̇ = (T(t) − T₀) · k_B / Q
```

where `Q = N_f · k_B · T₀ · τ²` is the thermostat mass.

In DMD, the thermostat is applied during the second half of the Velocity Verlet
step using a split-operator approach (half before half-kick, half after).

### 6.3 PME (Particle Mesh Ewald)

Electrostatics is the computationally most expensive term. Ewald summation
splits the Coulomb 1/r potential into two rapidly converging parts:

```
1/r = erfc(αr)/r + erf(αr)/r

E_Coul = ½ Σ' q_i q_j · erfc(α r_ij) / r_ij    [real space, cutoff]
       + ½ Σ' q_i q_j · erf(α r_ij) / r_ij     [reciprocal space, FFT]
       - (k_e · α / √π) · Σ q_i²                [self-energy, subtracted]
```

**PME** implements the reciprocal part in three steps:

1. **Charge spreading**: point charges are interpolated onto a regular 3D grid
   using B-splines of order n (typically 4, cubic). Each charge contributes to
   a small grid region (n³ points).

2. **3D FFT**: the charge grid is transformed to Fourier space via 3D FFT
   (r2c half-complex), multiplied by the Green's function
   `G(k) = 4π · exp(−k²/4α²) / k² · B(k)` (where B(k) is the B-spline
   correction), and back-transformed (c2r).

3. **Back-interpolation**: the reciprocal grid potential is re-interpolated
   at atomic positions using B-spline derivatives, yielding the forces.

**Complexity**: O(N log N) via FFT, compared to O(N²) for direct summation.
For N > 2000, PME is faster than direct sum at equal accuracy.

**Exclusion correction**: PME automatically includes all pairs, including 1-2
and 1-3 pairs that must be excluded. The `CoulombExclusion` component subtracts
the `q_i·q_j/r` term for excluded pairs after the PME calculation.

### 6.4 Verlet List

For short-range interactions (LJ and direct Coulomb), DMD maintains a
**Verlet list** (neighbor list within cutoff + skin distance). The list is
rebuilt every 10 steps to balance correctness (atoms cannot travel more than
the skin distance in 10 steps) and computational cost (O(N²) to rebuild,
O(N·M) to compute forces).

### 6.5 Steepest Descent Minimization

```python
info = dmd.minimize(system_data, n_steps=500, step_size=0.001, energy_tol=0.01)
```

The algorithm:
1. Compute forces `F = −∇E(pos)`
2. Move atoms along the force: `x ← x + α · F/|F|`
3. If energy decreased (`E_new < E_old`): accept, double α
4. If energy increased: reject, halve α, retry
5. Repeat until convergence (`|ΔE| < tol`) or `n_steps`

---

## 7. The H5MD Format

H5MD (**HDF5 Molecular Dynamics**) is an open standard for storing MD trajectories
based on HDF5. It was developed by the community as an open alternative to
proprietary or legacy formats.

### 8.1 Why H5MD over DCD, XTC, or TRR?

| Feature | H5MD | DCD (CHARMM) | XTC (GROMACS) | TRR (GROMACS) |
|---|---|---|---|---|
| **Open** | Yes, documented standard | No, proprietary binary | No, GROMACS-specific | No, GROMACS-specific |
| **Self-describing** | Yes (HDF5 contains metadata) | No | No | No |
| **Partial access** | Yes (select positions, frames) | Yes (but complex) | No (decompresses) | Yes |
| **Multiple datatypes** | Positions, velocities, forces, box, parameters | Positions only | Positions only | Positions + velocities + forces |
| **Extensible** | Yes (chunked datasets) | Yes (append) | No (rewrites) | Yes (append) |
| **Interoperable** | MDAnalysis, VMD, h5py, Julia, MATLAB | CHARMM, MDAnalysis | GROMACS, MDAnalysis, VMD | GROMACS, MDAnalysis |
| **Compression** | Yes (HDF5 filters: gzip, szip, lzf) | No | Yes (proprietary) | No |

**Key advantages of H5MD:**
- **Portable**: an H5MD file can be read by Python (h5py), VMD, MDAnalysis,
  Julia, MATLAB, and any language with HDF5 bindings.
- **Self-describing**: dataset names and attributes (units, time steps) are
  embedded in the file itself. No separate parameter file needed.
- **Random access**: read only the last 100 frames from a 1-million-frame file
  without decompressing everything.
- **Single file**: trajectory, box, time, and step are in the same file.
  Energies can be in a separate file or the same one.

### 8.2 H5MD File Structure

```
trajectory.h5
├── /h5md/
│   └── version = [1, 1]
├── /particles/
│   └── /atoms/
│       ├── /position/
│       │   ├── value     [n_frames, n_atoms, 3]  double
│       │   ├── time      [n_frames]               double
│       │   ├── step      [n_frames]               int
│       │   └── specification (attributes: units = "nm", ...)
│       ├── /velocity/
│       │   ├── value     [n_frames, n_atoms, 3]  double
│       │   ├── time      [n_frames]               double
│       │   └── step      [n_frames]               int
│       └── /box/
│           └── /edges/
│               ├── value [n_frames, 3, 3]         double
│               ├── time  [n_frames]               double
│               └── step  [n_frames]               int
```

The `value` datasets are **chunked** and **extendible**: new frames are appended
with `H5Dset_extent` + hyperslab, without recreating the file.

### 8.3 Reading H5MD with Python

```python
import h5py
import numpy as np

with h5py.File("traj.h5", "r") as f:
    pos = f["/particles/atoms/position/value"][:]    # (n_frames, n_atoms, 3)
    step = f["/particles/atoms/position/step"][:]    # (n_frames,)
    time = f["/particles/atoms/position/time"][:]    # (n_frames,)
    box = f["/particles/atoms/box/edges/value"][:]   # (n_frames, 3, 3)

print(f"Frames: {pos.shape[0]}, Atoms: {pos.shape[1]}")
```

With MDAnalysis:
```python
import MDAnalysis as mda
u = mda.Universe("traj.h5")
for ts in u.trajectory:
    print(ts.frame, ts.time, ts.positions.shape)
```

---

## 8. Checkpoint and Restart

### 9.1 Python Checkpoint (JSON)

```python
import dmd
import json

# After simulation
dmd.write_checkpoint("checkpoint.json", engine, system_data)

# To resume
engine = dmd.run(checkpoint="checkpoint.json", config_data=config)
```

The JSON checkpoint contains:
- `format: "dmd_checkpoint"`, `version: 1`
- `system`: the entire `system.json` (topology + parameters), so restart
  **does not depend on the original force field** (self-contained)
- `state`: current step, time, positions, velocities, box_size, energy

### 9.2 C++ Binary Checkpoint

During simulation, if `checkpoint_interval > 0`, the C++ engine writes a
binary checkpoint (`magic = 0x444D4450`) with:
- Full state: step, time, PE, positions, velocities, forces, masses, charges,
  atom types
- Verlet list state: `step_since_rebuild` (guarantees **bit-exact** restart —
  the trajectory after restart is identical to continuous simulation)

### 9.3 Complete Checkpoint Workflow

```python
import dmd, json

# Step 1: fresh start
config = { ... }  # config.json dict
system = { ... }  # system.json dict
engine = dmd.run(system_data=system, config_data=config)

# Step 2: continue from checkpoint
engine2 = dmd.run(checkpoint="checkpoint.bin", config_data=config)

# Step 3: second continuation
engine3 = dmd.run(checkpoint="checkpoint.bin", config_data=config)
```

### 9.4 Minimization + Checkpoint + Production

```python
# 1. Minimize
info = dmd.minimize(system, n_steps=500, output={"checkpoint_path": "min.json"})

# 2. Equilibrate NVT
config_nvt = { ... thermostat: nose-hoover ... }
eng = dmd.run(checkpoint="min.json", config_data=config_nvt)

# 3. NPT production
config_npt = { ... thermostat: nose-hoover, barostat: andersen ... }
eng2 = dmd.run(checkpoint="checkpoint.bin", config_data=config_npt)
```

---

## 9. Appendix: Best Practices

### Choosing the Time Step (dt)

| System | Recommended dt | Notes |
|---|---|---|
| All-atom (no H constraints) | 0.001 ps (1 fs) | Safe for any system |
| With C–H, O–H, N–H constraints | 0.002 ps (2 fs) | Standard biomolecular MD |
| Coarse-grained | 0.005–0.020 ps | Depends on CG resolution |

### Choosing Cutoffs

| Interaction | Recommended cutoff | Notes |
|---|---|---|
| Lennard-Jones | 1.0–1.4 nm | 1.2 nm is the biomolecular standard |
| Coulomb (direct) | 1.0–1.4 nm | Must match LJ cutoff |
| PME grid spacing | 0.10–0.15 nm | Finer = more accurate, slower |
| PME spline order | 4 | Cubic = best quality/cost ratio |

### Coupling Constants

| Purpose | τ thermostat | τ barostat |
|---|---|---|
| Rapid equilibration | 0.1 ps | 1.0 ps |
| NVT production (Nose-Hoover) | 1.0–2.0 ps | — |
| NPT production | 1.0–2.0 ps | 5.0–10.0 ps |

### Performance

- **LJ + direct Coulomb**: O(N·M), dominant for small systems (< 2000 atoms)
- **PME**: O(N log N), faster than direct sum for N > 2000
- **Verlet list**: rebuild every 10 steps, O(N²) cost but amortized
- **ns/day**: for 1566 atoms (TIP3P, 2.5 nm box, NVT, PME): ~1.2 ns/day
  on one Apple Silicon M3 core

### Debugging and Validation

1. **Energy conservation**: in NVE, total energy (E_pot + E_kin) must oscillate
   without drift. Drift > 0.5% over 5000 steps indicates a problem.
2. **Bit-exact restart**: after a checkpoint, the trajectory must be identical to
   the continuous one. DMD guarantees this for the binary format.
3. **Non-NaN forces**: check that forces don't reach 10¹⁵ — a typical sign of
   wrong coordinates (wrong units, swapped columns).
4. **Temperature**: after NVT equilibration, temperature must oscillate around
   the target with fluctuations `√(2/N_f) · T₀` (~5% for 1000 atoms).

---

## 10. Appendix: Unit Conversion Summary

### CHARMM → DMD (performed automatically by charmm.py)

| Quantity | CHARMM (input) | DMD (internal) | Factor |
|---|---|---|---|
| σ (LJ) | Rmin/2 in Å | σ in nm | ÷10 × 2 ÷ 2^(1/6) |
| ε (LJ) | kcal/mol | kJ/mol | × 4.184 |
| k (bond) | kcal/(mol·Å²) | kJ/(mol·nm²) | × 4.184 × 100 |
| r₀ (bond) | Å | nm | ÷ 10 |
| k_θ (angle) | kcal/(mol·rad²) | kJ/(mol·rad²) | × 4.184 |
| θ₀ (angle) | degrees | radians | × π/180 |
| k_φ (dihedral) | kcal/mol | kJ/mol | × 4.184 |
| φ₀ (dihedral) | degrees | radians | × π/180 |

### PDB → DMD

PDB coordinates are in **Ångstroms** (traditionally). DMD expects **nm**.
The conversion (÷10) is the user's responsibility. The `from_pdb()` parser
keeps raw values.

### DMD → SI

| DMD | SI |
|---|---|
| 1 nm | 10⁻⁹ m |
| 1 kJ/mol | 1.66054 × 10⁻²¹ J per molecule |
| 1 ps | 10⁻¹² s |
| 1 kJ/(mol·nm) | 1.66054 × 10⁻¹² N per molecule |
