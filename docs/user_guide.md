# DMD User Guide

## 1. Introduzione

DMD (Delicious Molecular Dynamics) è un motore di simulazione di dinamica
molecolare ottimizzato per CPU. Permette di simulare sistemi atomici e
molecolari in diversi ensemble termodinamici (NVE, NVT, NPT).

Il software è composto da due livelli:

- **Core C++** — motore di simulazione vero e proprio: integratore Velocity
  Verlet, forze (Lennard-Jones, Coulomb, PME, bonded), termostati, barostati,
  checkpoint/restart.
- **Python API** — interfaccia utente principale: costruzione del sistema,
  conversione da formati standard (PDB, PSF, GRO), lancio della simulazione,
  accesso ai risultati.

Questa guida descrive il flusso di lavoro completo dalla preparazione del
sistema alla simulazione, usando esclusivamente il layer Python.

### 1.1 Ensemble supportati

| Ensemble | Descrizione |
|---|---|
| **NVE** | Microcanonico: energia, volume e numero di particelle costanti |
| **NVT** | Canonico: temperatura costante tramite termostato |
| **NPT** | Isobaro-isotermo: pressione e temperatura costanti |

---

## 2. Installazione

### 2.1 Dipendenze

| Dipendenza | Versione minima | Note |
|---|---|---|
| Python | 3.10 | 3.14 raccomandato |
| pip | — | — |
| CMake | 3.20 | Solo per compilazione del modulo C++ |
| Compilatore C++17 | gcc ≥ 9, clang ≥ 12 | Apple clang incluso in Xcode |
| FFTW | 3.3 | Opzionale, solo per PME |

### 2.2 Build del modulo C++

```bash
git clone <url-del-repository>
cd dmd

# Crea directory di build e configura
cmake -B build -DCMAKE_BUILD_TYPE=Release

# Compila tutto (C++ tests + modulo Python)
cmake --build build -j8

# (Opzionale) Esegui i test C++
ctest --test-dir build
```

Il modulo Python `_dmd_core.cpython-XXX-darwin.so` viene automaticamente
piazzato in `python/dmd/` e importato dal package Python.

### 2.3 Verifica installazione

```bash
cd python
PYTHONPATH=$(pwd) python3 -c "import dmd; print('OK')"
```

Se il comando stampa `OK` senza errori, l'installazione è completa.

---

## 3. Panoramica del flusso di lavoro

Una simulazione DMD si articola in 4 passi:

```
┌─────────────┐    ┌──────────────┐    ┌──────────────┐    ┌──────────┐
│ 1. Prepara   │ -> │ 2. Costruisci │ -> │ 3. Configura  │ -> │ 4. Eseguí │
│    sistema   │    │  system.json  │    │  config.json  │    │ dmd.run() │
│ (PDB/GRO/   │    │               │    │               │    │           │
│  PSF)       │    │               │    │               │    │           │
└─────────────┘    └──────────────┘    └──────────────┘    └──────────┘
```

1. **Prepara il sistema**: puoi usare un file PDB, GRO o PSF, oppure creare
   un file `system.json` a mano.
2. **Costruisci il sistema**: il `SystemBuilder` legge `system.json` e
   produce un oggetto `SimulationConfig` con topology, posizioni e parametri
   di force field.
3. **Configura la simulazione**: un file `config.json` specifica tutti i
   parametri di run (dt, n_steps, termostato, barostato, ecc.).
4. **Esegui**: `dmd.run()` costruisce il motore, esegue la simulazione e
   restituisce un oggetto `Engine` con i risultati.

---

## 4. Formato system.json

Il file `system.json` descrive il sistema atomico: scatola di simulazione,
posizioni, masse, cariche, tipi atomici, e parametri di force field.

### 4.1 Struttura generale

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
        "lennard_jones": {
            "pairs": [
                {"type_i": "Ar", "type_j": "Ar", "sigma": 0.3405, "epsilon": 0.997}
            ]
        },
        "bonds": [
            {"i": 0, "j": 1, "k": 1000.0, "r0": 0.15}
        ]
    }
}
```

### 4.2 Sezione `cell`

| Campo | Tipo | Obbligatorio | Descrizione |
|---|---|---|---|
| `positions` | array di [x,y,z] | sì | Coordinate atomiche in nm |
| `box_size` | float | sì | Lato della scatola cubica in nm |

### 4.3 Sezione `atoms`

| Campo | Tipo | Obbligatorio | Descrizione |
|---|---|---|---|
| `mass` | array di float | sì | Masse atomiche in g/mol |
| `charge` | array di float | no (default 0.0) | Cariche atomiche in e⁻ |
| `type` | array di string | sì | Nomi dei tipi atomici |

I tre array devono avere la stessa lunghezza (numero di atomi).

### 4.4 Sezione `force_field`

#### Lennard-Jones

```json
"lennard_jones": {
    "pairs": [
        {"type_i": "Ar", "type_j": "Ar", "sigma": 0.3405, "epsilon": 0.997},
        {"type_i": "OW", "type_j": "OW", "sigma": 0.3166, "epsilon": 0.650},
        {"type_i": "Ar", "type_j": "OW", "sigma": 0.3285, "epsilon": 0.805}
    ]
}
```

I parametri sigma/epsilon vengono espansi in matrici `n_types × n_types`.
`n_types` è il numero di tipi atomici distinti nella sezione `atoms.type`.

Se omesso, Lennard-Jones non viene calcolato.

#### Legami (bonds)

```json
"bonds": [
    {"i": 0, "j": 1, "k": 1000.0, "r0": 0.15},
    {"i": 1, "j": 2, "k": 1000.0, "r0": 0.15}
]
```

| Campo | Tipo | Obbligatorio | Descrizione |
|---|---|---|---|
| `i`, `j` | int | sì | Indici atomici (0-based) |
| `k` | float | sì | Costante di forza in kJ/(mol·nm²) |
| `r0` | float | sì | Lunghezza di equilibrio in nm |

---

## 5. Formato config.json

Il file `config.json` specifica tutti i parametri della simulazione. Usa
uno schema **stretto**: ogni chiave deve essere presente; chiavi sconosciute
producono errore. Questo garantisce che non ci siano valori di default
nascosti.

Puoi generare un template completo con:

```python
import dmd
print(dmd.generate_config_template())
```

### 5.1 Schema completo

```json
{
    "run": {
        "dt": 0.002,
        "n_steps": 10000,
        "init_temperature": 300.0,
        "gen_vel": false,
        "seed": 42
    },
    "output": {
        "trajectory_path": "",
        "nstxout": 100,
        "nstvout": 100,
        "nstenergy": 10,
        "energy_path": "",
        "checkpoint_interval": 0,
        "checkpoint_path": ""
    },
    "lj": {
        "cutoff": 1.2
    },
    "electrostatics": {
        "coulomb_type": "cutoff",
        "cutoff": 1.2,
        "pme_order": 4,
        "pme_grid_spacing": 0.12,
        "ewald_coeff": 3.5
    },
    "thermostat": {
        "type": "none",
        "temperature": 300.0,
        "tau": 0.1,
        "frequency": 10.0
    },
    "barostat": {
        "type": "none",
        "pressure": 1.0,
        "tau": 1.0,
        "compressibility": 4.5e-5
    },
    "constraints": {
        "type": "none",
        "tolerance": 1e-6
    }
}
```

### 5.2 Reference sezione `run`

| Chiave | Tipo | Descrizione |
|---|---|---|
| `dt` | float | Passo di integrazione in picosecondi |
| `n_steps` | int | Numero totale di passi |
| `init_temperature` | float | Temperatura di generazione velocità iniziali (K) |
| `gen_vel` | bool | Se true, genera velocità Maxwell-Boltzmann all'avvio |
| `seed` | int | Seme del generatore casuale per velocità iniziali |

### 5.3 Reference sezione `output`

| Chiave | Tipo | Descrizione |
|---|---|---|
| `trajectory_path` | string | Percorso file H5MD della traiettoria (vuoto = nessun output) |
| `nstxout` | int | Scrivi posizioni ogni N passi |
| `nstvout` | int | Scrivi velocità ogni N passi (0 = salta) |
| `nstenergy` | int | Scrivi energie ogni N passi |
| `energy_path` | string | Percorso file H5MD delle osservabili |
| `checkpoint_interval` | int | Checkpoint ogni N passi (0 = disabilitato) |
| `checkpoint_path` | string | Percorso file checkpoint |

### 5.4 Reference sezione `lj`

| Chiave | Tipo | Descrizione |
|---|---|---|
| `cutoff` | float | Raggio di cutoff Lennard-Jones in nm |

### 5.5 Reference sezione `electrostatics`

| Chiave | Tipo | Descrizione |
|---|---|---|
| `coulomb_type` | string | Metodo: `"cutoff"`, `"pme"`, o `"direct"` |
| `cutoff` | float | Raggio di cutoff Coulomb in nm |
| `pme_order` | int | Ordine B-spline per PME |
| `pme_grid_spacing` | float | Spaziatura griglia PME in nm |
| `ewald_coeff` | float | Coefficiente Ewald α in nm⁻¹ |

### 5.6 Reference sezione `thermostat`

| Chiave | Tipo | Descrizione |
|---|---|---|
| `type` | string | `"none"`, `"berendsen"`, `"nose-hoover"`, `"andersen"` |
| `temperature` | float | Temperatura target (K) |
| `tau` | float | Costante di accoppiamento (ps) — Berendsen e Nose-Hoover |
| `frequency` | float | Frequenza di collisione (ps⁻¹) — Andersen |

### 5.7 Reference sezione `barostat`

| Chiave | Tipo | Descrizione |
|---|---|---|
| `type` | string | `"none"`, `"berendsen"`, `"andersen"` |
| `pressure` | float | Pressione target (bar) |
| `tau` | float | Costante di accoppiamento (ps) |
| `compressibility` | float | Comprimibilità isoterma (bar⁻¹) — Berendsen |

### 5.8 Reference sezione `constraints`

| Chiave | Tipo | Descrizione |
|---|---|---|
| `type` | string | `"none"` (unica opzione supportata) |
| `tolerance` | float | Tolleranza del costruttore di vincoli |

---

## 6. Tutorial: Argon NVE

Questo tutorial mostra una simulazione NVE completa di 32 atomi di Argon
in una scatola cubica da 1.5 nm, partendo da zero.

### 6.1 Crea system.json

Crea un file `system.json` con il seguente contenuto:

```json
{
    "cell": {
        "positions": [
            [0.000, 0.000, 0.000],
            [0.375, 0.375, 0.000],
            [0.375, 0.000, 0.375],
            [0.000, 0.375, 0.375],
            [0.750, 0.750, 0.000],
            [1.125, 1.125, 0.000],
            [1.125, 0.750, 0.375],
            [0.750, 1.125, 0.375],
            [0.750, 0.000, 0.750],
            [1.125, 0.375, 0.750],
            [1.125, 0.000, 1.125],
            [0.750, 0.375, 1.125],
            [0.000, 0.750, 0.750],
            [0.375, 1.125, 0.750],
            [0.375, 0.750, 1.125],
            [0.000, 1.125, 1.125],
            [0.000, 0.000, 0.750],
            [0.375, 0.375, 0.750],
            [0.375, 0.000, 1.125],
            [0.000, 0.375, 1.125],
            [0.750, 0.750, 0.750],
            [1.125, 1.125, 0.750],
            [1.125, 0.750, 1.125],
            [0.750, 1.125, 1.125],
            [0.750, 0.000, 0.000],
            [1.125, 0.375, 0.000],
            [1.125, 0.000, 0.375],
            [0.750, 0.375, 0.375],
            [0.000, 0.750, 0.000],
            [0.375, 1.125, 0.000],
            [0.375, 0.750, 0.375],
            [0.000, 1.125, 0.375]
        ],
        "box_size": 1.5
    },
    "atoms": {
        "mass": [39.948, 39.948, 39.948, 39.948, 39.948, 39.948, 39.948, 39.948, 39.948, 39.948, 39.948, 39.948, 39.948, 39.948, 39.948, 39.948, 39.948, 39.948, 39.948, 39.948, 39.948, 39.948, 39.948, 39.948, 39.948, 39.948, 39.948, 39.948, 39.948, 39.948, 39.948, 39.948],
        "charge": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        "type": ["Ar", "Ar", "Ar", "Ar", "Ar", "Ar", "Ar", "Ar", "Ar", "Ar", "Ar", "Ar", "Ar", "Ar", "Ar", "Ar", "Ar", "Ar", "Ar", "Ar", "Ar", "Ar", "Ar", "Ar", "Ar", "Ar", "Ar", "Ar", "Ar", "Ar", "Ar", "Ar"]
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

> Le coordinate rappresentano un reticolo FCC 2×2×2 con parametro di reticolo
> 0.75 nm (32 atomi in 1.5 nm³).

### 6.2 Crea config.json

Crea un file `config.json` per una simulazione NVE a 85 K:

```json
{
    "run": {
        "dt": 0.002,
        "n_steps": 5000,
        "init_temperature": 85.0,
        "gen_vel": true,
        "seed": 42
    },
    "output": {
        "trajectory_path": "",
        "nstxout": 0,
        "nstvout": 0,
        "nstenergy": 0,
        "energy_path": "",
        "checkpoint_interval": 0,
        "checkpoint_path": ""
    },
    "lj": {
        "cutoff": 1.0
    },
    "electrostatics": {
        "coulomb_type": "cutoff",
        "cutoff": 1.0,
        "pme_order": 4,
        "pme_grid_spacing": 0.12,
        "ewald_coeff": 3.5
    },
    "thermostat": {
        "type": "none",
        "temperature": 85.0,
        "tau": 0.1,
        "frequency": 10.0
    },
    "barostat": {
        "type": "none",
        "pressure": 1.0,
        "tau": 1.0,
        "compressibility": 4.5e-5
    },
    "constraints": {
        "type": "none",
        "tolerance": 1e-6
    }
}
```

### 6.3 Lancia la simulazione

```python
import dmd
import json

with open("system.json") as f:
    system_data = json.load(f)

with open("config.json") as f:
    config_data = json.load(f)

engine = dmd.run(system_data=system_data, config_data=config_data)

print(f"Passo finale: {engine.current_step}")
print(f"Tempo finale: {engine.current_time:.4f} ps")
print(f"Energia potenziale: {engine.potential_energy:.4f} kJ/mol")
print(f"Posizioni finali:\n{engine.positions}")
```

Output atteso (valori indicativi):

```
Passo finale: 5000
Tempo finale: 10.0000 ps
Energia potenziale: -85.2341 kJ/mol
Posizioni finali:
[[...]]
```

### 6.4 Analisi rapida

```python
import numpy as np

pos = np.array(engine.positions)
vel = np.array(engine.velocities)

# Energia cinetica
masses = np.array(system_data["atoms"]["mass"])
ekin = 0.5 * np.sum(masses * np.sum(vel**2, axis=1))
print(f"Energia cinetica: {ekin:.4f} kJ/mol")

# Temperatura istantanea
N = len(masses)
kb = 0.008314462618  # kJ/(mol·K)
T = 2 * ekin / (3 * N * kb)
print(f"Temperatura: {T:.2f} K")

# Raggio di girazione
com = np.sum(masses[:, None] * pos, axis=0) / np.sum(masses)
rg2 = np.sum(masses * np.sum((pos - com)**2, axis=1)) / np.sum(masses)
print(f"Raggio di girazione^2: {rg2:.4f} nm^2")
```

---

## 7. Convertitori di formato

Se hai un sistema in formato PDB, GRO o PSF, DMD può convertirlo
automaticamente in system.json.

### 7.1 Da PDB

```python
import dmd

system = dmd.from_pdb("input.pdb")

# system è un dict compatibile con SystemBuilder
builder = dmd.SystemBuilder(system)
cfg = builder.build()
```

Il parser PDB legge:
- `CRYST1` → `cell.box_size`
- `ATOM` / `HETATM` → posizioni e tipo atomico (dal campo elemento, colonne 76-78)
- La massa è stimata dal tipo atomico

### 7.2 Da GRO

```python
import dmd

system = dmd.from_gro("input.gro")
builder = dmd.SystemBuilder(system)
```

Il parser GRO legge:
- Prima riga: commento (ignorato)
- Seconda riga: numero atomi
- Righe atomi: residuo, nome atomo, posizioni x/y/z
- Ultima riga: dimensioni scatola

### 7.3 Da PSF

```python
import dmd
from dmd.convert.psf import from_psf

psf_data = from_psf("structure.psf")
# Aggiunge bonds/angles/dihedrals a un system.json esistente
with open("system.json") as f:
    system = json.load(f)
system["force_field"].update(psf_data["force_field"])
```

Il parser PSF legge le sezioni:
- `!NBOND` → legami
- `!NTHETA` → angoli
- `!NPHI` → diedri
- `!NIMPHI` → impropri

### 7.4 Workflow combinato

```python
import dmd
import json

# Leggi coordinate da PDB
system = dmd.from_pdb("protein.pdb")

# Aggiungi topology da PSF
psf = dmd.from_psf("protein.psf")
system["force_field"].update(psf["force_field"])

# Configura e lancia
config = {
    "run": {"dt": 0.002, "n_steps": 10000, "init_temperature": 300.0,
            "gen_vel": True, "seed": 42},
    "output": {"trajectory_path": "", "nstxout": 0, "nstvout": 0,
               "nstenergy": 0, "energy_path": "",
               "checkpoint_interval": 0, "checkpoint_path": ""},
    "lj": {"cutoff": 1.2},
    "electrostatics": {"coulomb_type": "cutoff", "cutoff": 1.2,
                       "pme_order": 4, "pme_grid_spacing": 0.12,
                       "ewald_coeff": 3.5},
    "thermostat": {"type": "none", "temperature": 300.0, "tau": 0.1,
                   "frequency": 10.0},
    "barostat": {"type": "none", "pressure": 1.0, "tau": 1.0,
                 "compressibility": 4.5e-5},
    "constraints": {"type": "none", "tolerance": 1e-6},
}

engine = dmd.run(system_data=system, config_data=config)
print(f"Energia finale: {engine.potential_energy:.4f} kJ/mol")
```

---

## 8. Riferimento API Python

### 8.1 `dmd.SystemBuilder`

Costruisce un `SimulationConfig` a partire da un dict system.json.

```python
from dmd import SystemBuilder

# Da dict
builder = SystemBuilder(system_dict)
cfg = builder.build()

# Da stringa JSON
builder = SystemBuilder(json_string)

# Concatena config.json
builder.apply_config_json(config_dict)
builder.apply_config_json("config.json")
```

### 8.2 `dmd.run()`

Funzione orchestratore: accetta system.json e config.json in varie forme
e restituisce un motore già eseguito.

```python
engine = dmd.run(
    system_data=system_dict,          # dict system.json
    config_data=config_dict,          # dict config.json
)
```

### 8.3 `dmd.SimulationConfig`

Configurazione completa della simulazione. Accessibile come oggetto Python
con proprietà numpy array.

```python
cfg = dmd.SimulationConfig()

# Scalari
cfg.dt = 0.002
cfg.n_steps = 10000
cfg.box_size = 3.0
cfg.lj_cutoff = 1.2
cfg.init_temperature = 300.0
cfg.gen_vel = True
cfg.seed = 42

# Array (numpy)
import numpy as np
cfg.mass = np.array([39.948, 39.948])
cfg.charges = np.array([0.0, 0.0])
cfg.pos_x = np.array([0.0, 1.0])
cfg.pos_y = np.array([0.0, 0.0])
cfg.pos_z = np.array([0.0, 0.0])
cfg.atom_types = np.array([0, 0], dtype=np.int32)

# Parametri LJ
cfg.lj_sigma = np.array([0.3405])        # n_types² elementi
cfg.lj_epsilon = np.array([0.997])        # n_types² elementi

# Legami
cfg.bond_i = np.array([0, 1], dtype=np.int32)
cfg.bond_j = np.array([1, 2], dtype=np.int32)
cfg.bond_k = np.array([1000.0, 1000.0])
cfg.bond_r0 = np.array([0.15, 0.15])
```

### 8.4 `dmd.Engine`

Motore di simulazione eseguito. Fornisce accesso ai risultati.

```python
engine = dmd.Engine(cfg)
engine.run()

# Proprietà in lettura
engine.current_step      # int — passo corrente
engine.current_time      # float — tempo corrente (ps)
engine.positions         # ndarray (n_atoms, 3) — posizioni
engine.velocities        # ndarray (n_atoms, 3) — velocità
engine.potential_energy  # float — energia potenziale (kJ/mol)
```

### 8.5 Funzioni libere

```python
import dmd

# Template config
template = dmd.generate_config_template()
```

### 8.6 Convertitori

```python
from dmd import from_pdb, from_psf, from_gro

system = from_pdb("file.pdb")       # -> dict system.json
system = from_gro("file.gro")       # -> dict system.json
psf_data = from_psf("file.psf")     # -> dict {force_field: {...}}
```

---

## 9. Simulazioni avanzate

### 9.1 NVT con termostato Berendsen

```python
config = {
    "run": {"dt": 0.002, "n_steps": 10000, "init_temperature": 300.0,
            "gen_vel": True, "seed": 42},
    # ... altri campi ...
    "thermostat": {
        "type": "berendsen",
        "temperature": 300.0,
        "tau": 0.1
    },
    "barostat": {"type": "none", "pressure": 1.0, "tau": 1.0,
                 "compressibility": 4.5e-5},
    # ...
}
```

### 9.2 NPT con Berendsen (termostato + barostato)

```python
config = {
    "run": {"dt": 0.002, "n_steps": 20000, "init_temperature": 120.0,
            "gen_vel": True, "seed": 42},
    # ... altri campi ...
    "thermostat": {
        "type": "berendsen",
        "temperature": 120.0,
        "tau": 0.1
    },
    "barostat": {
        "type": "berendsen",
        "pressure": 100.0,
        "tau": 1.0,
        "compressibility": 4.5e-5
    },
    # ...
}
```

### 9.3 Con PME ( Particle Mesh Ewald )

```python
config = {
    # ... run, output, lj ...
    "electrostatics": {
        "coulomb_type": "pme",
        "cutoff": 1.2,
        "pme_order": 4,
        "pme_grid_spacing": 0.12,
        "ewald_coeff": 3.5
    },
    # ...
}
```

### 9.4 Checkpoint e restart

```python
config = {
    # ... run ...
    "output": {
        # ...
        "checkpoint_interval": 500,
        "checkpoint_path": "checkpoint.bin"
    },
    # ...
}
```

Per riprendere da un checkpoint, usa `dmd.run(checkpoint="checkpoint.json", config_data=config)`

> Il checkpoint C++ è bit-exact: riprendendo da un checkpoint si ottiene
> la stessa traiettoria di una simulazione continua. Le variabili estese
> (Nose-Hoover ζ, Andersen piston) non sono ancora checkpointate.

---

## 10. Appendice A: Ensemble — teoria

### 10.1 Velocity Verlet

DMD usa l'integratore Velocity Verlet, un algoritmo simplettico (conserva
l'energia) e reversibile nel tempo.

```
v(t + ½Δt) = v(t) + (Δt/2) · a(t)          # half kick
r(t + Δt)  = r(t) + Δt · v(t + ½Δt)        # advance
a(t + Δt)  = F(t + Δt) / m                 # forze
v(t + Δt)  = v(t + ½Δt) + (Δt/2) · a(t + Δt)  # half kick
```

### 10.2 Termostato Berendsen

Scala le velocità per guidare la temperatura istantanea verso il target:

```
λ = √(1 + Δt/τ · (T0/T(t) - 1))
v_i ← λ · v_i
```

- `τ` piccolo → accoppiamento forte (utile per equilibrazione)
- Non produce un ensemble canonico esatto

### 10.3 Termostato Nose-Hoover

Estende il sistema con un coefficiente d'attrito ζ:

```
ζ ← ζ + Δt/2 · (T(t) - T0) · k_B / Q
v_i ← v_i · exp(-ζ · Δt/2)
Q = N_f · k_B · T0 · τ²
```

Produce una distribuzione canonica (NVT) esatta.

### 10.4 Termostato Andersen

Sostituisce casualmente le velocità delle particelle con estrazioni
Maxwell-Boltzmann:

```
P(collisione) = ν · Δt
v_i ← √(k_B · T0 / m_i) · N(0, 1)   per ogni componente
```

Produce un ensemble canonico ma interrompe la dinamica (traiettorie non
continue).

### 10.5 Barostato Berendsen

Scala il volume della scatola per guidare la pressione verso il target:

```
μ = 1 - Δt/τ_P · β · (P0 - P(t))
r_i ← μ^(1/3) · r_i
box ← μ^(1/3) · box
```

### 10.6 Barostato Andersen

Usa un grado di libertà aggiuntivo ( pistone ϵ ) con lagrangiana estesa:

```
ϵ̇ ← ϵ̇ + Δt/2 · 3V · (P(t) - P0) / W
r_i ← r_i · exp(ϵ̇ · Δt)
box ← box · exp(ϵ̇ · Δt)
W = N_f · k_B · P0 · τ²
```

---

## 11. Appendice B: Forze

### 11.1 Lennard-Jones

```
E_LJ(rij) = 4ε_ij [ (σ_ij/rij)¹² - (σ_ij/rij)⁶ ]    per r < r_cutoff
F_ij = 24 ε_ij (2 (σ_ij/rij)¹² - (σ_ij/rij)⁶) / rij² · r_ij
```

Per sistemi multi-tipo, i parametri incrociati sono precalcolati in matrici
`n_types × n_types`. La lista di vicini di Verlet viene ricostruita ogni
10 passi.

### 11.2 Coulomb (direct sum)

```
E_Coul(rij) = k_e · qi · qj · erfc(α rij) / rij    per r < r_cutoff
```

k_e = 138.935456 kJ·mol⁻¹·nm·e⁻², α è il coefficiente Ewald.

### 11.3 Coulomb (PME)

L'energia elettrostatica totale è divisa in tre parti:

```
E_Coul = E_real + E_recip + E_self
```

- **Reale**: somma diretta con erfc screening (come sopra, cutoff)
- **Reciproca**: cariche interpolate su griglia 3D con B-spline, FFT 3D,
  funzione di Green in spazio di Fourier, back-interpolazione forze
- **Self-energy**: correzione −k_e · α / √π · Σ q_i²

### 11.4 Bonded

**Legame armonico:** E_bond = ½ k_b (r - r0)²

r = |ri - rj|, k_b costante di forza, r0 lunghezza equilibrio.

**Angolo armonico:** E_angle = ½ k_θ (θ - θ0)²

θ = angolo tra atomi i-j-k.

**Diedro periodico:** E_dihedral = k_φ (1 + cos(n·φ - φ0))

φ = diedro tra atomi i-j-k-l, n periodicità, φ0 fase.

---

## 12. Appendice D: Best practices

### 12.1 Scelta del passo temporale

| Sistema | dt raccomandato |
|---|---|
| Tutti gli atomi (standard) | 0.002 ps (2 fs) |
| Con atomi H senza constraints | 0.001 ps |
| Coarse-grained | 0.005 ps o più |

### 12.2 Scelta dei cutoff

| Interazione | Cutoff raccomandato |
|---|---|
| Lennard-Jones | 1.0 – 1.4 nm |
| Coulomb (cutoff) | Come LJ |
| PME grid spacing | 0.10 – 0.15 nm, spline order 4 |

### 12.3 Costanti di accoppiamento termostato/barostato

| Scopo | τ termostato | τ barostato |
|---|---|---|
| Equilibrazione forte | 0.1 ps | 1.0 ps |
| Produzione debole | 1.0 ps | 5.0 ps |

Per produzione NPT: Nose-Hoover + Andersen barostat (ensemble più accurato).

### 12.4 Checkpoint

- Scrivi checkpoint ogni 500 – 5000 passi
- Verifica sempre che restart riproduca traiettorie identiche

### 12.5 Performance

- LJ e Coulomb direct dominano il costo computazionale (loop a coppie)
- PME scala O(N log N) via FFT, più efficiente del direct sum per N > 2000
- Compila con `-O2` o `-O3` per produzione

---

## 13. Appendice E: Struttura del codice

```
dmd/
├── python/                    # Package Python
│   ├── dmd/
│   │   ├── __init__.py        # Esporta API pubblica
│   │   ├── core.cpp           # Bindings pybind11 → _dmd_core.so
│   │   ├── system.py          # SystemBuilder
│   │   ├── runner.py          # run() orchestrator
│   │   ├── config.py          # Config Parser (validate, apply)
│   │   └── convert/
│   │       ├── __init__.py
│   │       ├── pdb.py         # Parse PDB
│   │       ├── psf.py         # Parse PSF
│   │       └── gro.py         # Parse GRO
│   └── pyproject.toml
├── src/                       # Core C++
│   ├── sim/                   # SimulationEngine, checkpoin
│   ├── force/                 # Forze (LJ, Coulomb, PME, bonded)
│   ├── cell/                  # Cella di simulazione, neighbor list
│   ├── thermo/                # Termostati
│   └── baro/                  # Barostati
└── tests/                     # Test C++ (Google Test)
```
