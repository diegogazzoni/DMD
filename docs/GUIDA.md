# DMD — Guida Utente Completa

## 1. Introduzione alla Dinamica Molecolare

### 1.1 Cos'è la Dinamica Molecolare?

La **Dinamica Molecolare** (MD) è una tecnica di simulazione al computer che studia il
movimento di atomi e molecole nel tempo. Ogni atomo viene trattato come una particella
puntiforme con una massa e una carica; il suo moto è governato dalle leggi della meccanica
classica (seconda legge di Newton):

```
F = m · a
```

Integrando numericamente questa equazione per ogni atomo a ogni passo temporale
(Δt ~ 1–2 femtosecondi), si ottiene una **traiettoria**: l'evoluzione del sistema
nell'arco di nanosecondi o microsecondi.

La MD permette di calcolare proprietà termodinamiche (temperatura, pressione, energia),
strutturali (distribuzioni radiali, angoli di legame) e dinamiche (coefficienti di
diffusione, spettri vibrazionali) di sistemi che vanno da poche centinaia a milioni
di atomi.

### 1.2 Potenziali e Force Field

Le forze tra atomi derivano da **campi di forza** (force field), funzioni potenziali
semi-empiriche che approssimano l'energia potenziale del sistema come somma di
contributi:

```
E_tot = E_bond + E_angle + E_dihedral + E_LJ + E_Coulomb
```

| Contributo | Descrizione | Modello |
|---|---|---|
| **E_bond** | Legame chimico tra due atomi | Armonico: ½k(r − r₀)² |
| **E_angle** | Angolo tra tre atomi | Armonico: ½k_θ(θ − θ₀)² |
| **E_dihedral** | Rotazione attorno a un legame | Periodico: k_φ[1+cos(nφ−φ₀)] |
| **E_LJ** (Lennard-Jones) | Forze di van der Waals | 4ε[(σ/r)¹² − (σ/r)⁶] |
| **E_Coulomb** | Interazioni elettrostatiche | k_e · q_i · q_j / r |

Di questi, **bond/angle/dihedral** sono detti **forze di legame** (bonded) perché
coinvolgono atomi connessi da legami covalenti. **LJ e Coulomb** sono dette
**forze non di legame** (non-bonded) perché agiscono tra tutte le coppie di atomi,
indipendentemente dalla connettività.

### 1.3 Short-range vs Long-range

La distinzione più importante dal punto di vista computazionale è tra forze
a **corto raggio** e a **lungo raggio**:

- **Short-range** (LJ, Coulomb direct): decadono rapidamente con r. Si applica
  un **cutoff** (tipicamente 1.0–1.4 nm): oltre quella distanza il contributo
  è trascurabile e non viene calcolato. La complessità è O(N · M) dove M è il
  numero medio di vicini entro il cutoff (~50–500 atomi).

- **Long-range** (Coulomb): l'interazione elettrostatica decade come 1/r,
  troppo lentamente per essere troncata. Usare un cutoff causerebbe artefatti
  gravi. Si usano metodi come **Ewald summation** o **Particle Mesh Ewald (PME)**
  che separano il potenziale in una parte reale (calcolata in spazio diretto con
  cutoff) e una parte reciproca (calcolata in spazio di Fourier via FFT 3D).
  Complessità: O(N log N).

### 1.4 Vincoli (Constraints)

In MD classica, i legami che coinvolgono idrogeni (C–H, O–H, N–H) vibrano a
frequenze molto elevate (~3000 cm⁻¹). Per integrarle correttamente servirebbe
Δt < 0.5 fs, rallentando la simulazione. La soluzione è **vincolare** queste
lunghezze a un valore fisso usando algoritmi come **SHAKE** o **LINCS**, che
permettono di usare Δt = 2 fs anche in presenza di idrogeni.

I vincoli vengono applicati _dopo_ l'integrazione correggendo posizioni e/o velocità.

### 1.5 Ensemble termodinamici

| Ensemble | Costanti | Descrizione |
|---|---|---|
| **NVE** | Numero atomi, Volume, Energia | Microcanonico: sistema isolato, energia totale conservata |
| **NVT** | Numero atomi, Volume, Temperatura | Canonico: temperatura costante via termostato |
| **NPT** | Numero atomi, Pressione, Temperatura | Isobaro-isotermo: pressione e temperatura costanti |

In NVT si usa un **termostato** che scambia calore con un bagno termico virtuale.
In NPT si aggiunge un **barostato** che scala il volume per mantenere la pressione
target.

---

## 2. Architettura del Software

DMD è organizzato su due livelli:

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
│     Checkpoint (binario)                │  sim/checkpoint.cpp
├─────────────────────────────────────────┤
│   Cell / SystemData / Constants         │  core/*.cpp
├─────────────────────────────────────────┤
│            C++ Core Layer               │
└─────────────────────────────────────────┘
```

### Principi architetturali

- **Separazione netta**: C++ gestisce il loop critico (forze, integrazione, I/O
  traiettoria). Python gestisce setup, conversione formati, orchestrazione.
- **SystemData SoA**: tutte le proprietà atomiche sono in vettori paralleli
  (Structure of Arrays): `pos_x`, `pos_y`, `pos_z` invece di `pos[i].x`. Massimizza
  la località della cache e facilita l'uso di span.
- **ForceComponent**: ogni termine energetico è una classe che implementa
  `compute(span pos_x/y/z, cell, span forces_x/y/z, energy&)`. Si registrano
  all'avvio in un ordine fisso.
- **Config JSON strict**: il file config.json segue uno schema rigido: ogni
  chiave deve essere presente, chiavi sconosciute producono errore. Niente
  default nascosti.
- **H5MD nativo**: la traiettoria è in formato H5MD (HDF5), uno standard aperto.
  Nessun formato proprietario.

### Unità e conversioni

Tutte le quantità interne usano il sistema **DMD** (coerente con GROMACS):
- Lunghezza: **nm**
- Energia: **kJ/mol**
- Forza: **kJ/(mol·nm)**
- Tempo: **ps**
- Carica: **carica elementare e⁻**

I parser dei force field (es. CHARMM, che usa kcal/mol e Å) convertono
automaticamente all'import.

---

## 3. Formato system.json

Il file **system.json** descrive il sistema atomico: la scatola di simulazione,
le posizioni atomiche, le masse, le cariche, i tipi atomici e i parametri di
force field. Può essere creato a mano o generato dai convertitori (PDB, PSF, GRO).

### 3.1 Struttura generale

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

### 3.2 Sezione `cell`

| Campo | Tipo | Obbligatorio | Descrizione |
|---|---|---|---|
| `positions` | array di [x,y,z] | sì | Coordinate atomiche in nm |
| `box_size` | float | sì | Lato della scatola cubica in nm |

La scatola è attualmente **cubica** (un solo lato). In futuro celle tricline,
ortorombiche, dodecaedro rombico.

### 3.3 Sezione `atoms`

| Campo | Tipo | Obbligatorio | Descrizione |
|---|---|---|---|
| `mass` | array di float | sì | Masse atomiche in g/mol |
| `charge` | array di float | sì | Cariche atomiche in e⁻ |
| `type` | array di string | sì | Nomi dei tipi atomici (es. "OH2", "HW1", "Ar") |

I tre array devono avere la stessa lunghezza (N atomi). I nomi dei tipi
devono corrispondere a quelli usati in `force_field.lennard_jones.pairs`.

### 3.4 Sezione `force_field`

#### Lennard-Jones

```json
"lennard_jones": {
    "pairs": [
        {"type_i": "Ar", "type_j": "Ar", "sigma": 0.3405, "epsilon": 0.997},
        {"type_i": "OH2", "type_j": "OH2", "sigma": 0.3150, "epsilon": 0.636}
    ]
}
```

| Campo | Tipo | Obbligatorio | Descrizione |
|---|---|---|---|
| `type_i`, `type_j` | string | sì | Nomi dei tipi atomici |
| `sigma` | float | sì | Parametro LJ σ in **nm** |
| `epsilon` | float | sì | Parametro LJ ε in **kJ/mol** |

Attenzione: σ è il vero parametro LJ (posizione del minimo = 2^(1/6)·σ),
**non** Rmin/2 (il raggio di van der Waals). Se usate parametri CHARMM,
il parser `charmm.py` converte automaticamente:
- Rmin/2 (Å) → σ (nm): `σ = (2 · Rmin/2) / (10 · 2^(1/6))`
- ε (kcal/mol) → ε (kJ/mol): `ε × 4.184`

Le regole di combinazione per coppie miste (se non specificate esplicitamente)
seguono Lorentz-Berthelot:
- `σ_ij = (σ_i + σ_j) / 2`  (media aritmetica)
- `ε_ij = √(ε_i · ε_j)`     (media geometrica)

Le coppie vengono espanse in una matrice `n_types × n_types` all'avvio.

#### Legami (bonds)

```json
"bonds": [
    {"i": 0, "j": 1, "k": 502416.0, "r0": 0.09572}
]
```

| Campo | Tipo | Obbligatorio | Descrizione |
|---|---|---|---|
| `i`, `j` | int | sì | Indici atomici (0-based) |
| `k` | float | sì | Costante di forza in kJ/(mol·nm²) |
| `r0` | float | sì | Lunghezza di equilibrio in nm |

Energia: `E = ½ · k · (r − r₀)²` dove `r = |r_i − r_j|`.

#### Angoli (angles)

```json
"angles": [
    {"i": 0, "j": 1, "k": 2, "k_theta": 628.02, "theta0": 1.824}
]
```

| Campo | Tipo | Obbligatorio | Descrizione |
|---|---|---|---|
| `i`, `j`, `k` | int | sì | Indici atomici, j = atomo centrale |
| `k_theta` | float | sì | Costante di forza in kJ/(mol·rad²) |
| `theta0` | float | sì | Angolo di equilibrio in **radianti** |

Energia: `E = ½ · k_θ · (θ − θ₀)²`.

#### Diedri propri (dihedrals)

```json
"dihedrals": [
    {"i": 0, "j": 1, "k": 2, "l": 3, "k_phi": 5.0, "periodicity": 3, "phi0": 0.0}
]
```

| Campo | Tipo | Obbligatorio | Descrizione |
|---|---|---|---|
| `i`, `j`, `k`, `l` | int | sì | Indici atomici |
| `k_phi` | float | sì | Costante di forza in kJ/mol |
| `periodicity` | int | sì | Molteplicità n |
| `phi0` | float | sì | Fase in **radianti** |

Energia: `E = k_φ · [1 + cos(n · φ − φ₀)]`.

Multi-term diedri: più entry con la stessa quadrupla (i,j,k,l) ma diversi
`k_phi`, `periodicity`, `phi0`. Il parser CHARMM le gestisce automaticamente.

#### Improper (diedri impropri)

```json
"impropers": [
    {"i": 0, "j": 1, "k": 2, "l": 3, "k_phi": 10.0, "phi0": 3.14159}
]
```

Stessi campi dei diedri propri. Descrivono torsioni fuori dal piano
(es. atomo centrale planare in un gruppo carbonilico). Vengono calcolati
dallo stesso componente `PeriodicDihedral`.

### 3.5 Sezione `exclusions`

```json
"exclusions": [[0, 1], [0, 2], [1, 3]]
```

Coppie di atomi per cui **non** calcolare le interazioni non-bonded
(LJ + Coulomb). Tipicamente:
- **1-2 exclusion**: atomi legati direttamente
- **1-3 exclusion**: atomi separati da due legami (angolo)
- **1-4 exclusion** (opzionale): scala le interazioni diedro

Il parser PSF (`from_psf()`) genera automaticamente le esclusioni dalla
connettività (bond → 1-2, angle → 1-3) e può includere lo scaling 1-4
con `scale_14`. I parser PDB e GRO non generano esclusioni.

### 3.6 Generazione del system.json

#### Da PDB + PSF + Force Field (flusso completo)

```python
import dmd
import json

# 1. Leggi coordinate
system = dmd.from_pdb("struttura.pdb")

# 2. Leggi topology (legami, angoli, diedri, esclusioni)
psf = dmd.from_psf("struttura.psf", scale_14=1.0)
system["atoms"]["type"] = psf["atom_types"]
system["force_field"].update(psf["force_field"])

# 3. Carica parametri force field
ff_params = dmd.load_forcefield("parametri.prm", format="charmm")

# 4. Merge: abbina parametri ai tipi atomici del PSF
merged = dmd.merge_ff(system, ff_params)
system["force_field"].update(merged["force_field"])
system["exclusions"] = psf["exclusions"]

# 5. Salva
with open("system.json", "w") as f:
    json.dump(system, f, indent=2)
```

#### Da GRO (solo coordinate, nessun bonded)

```python
system = dmd.from_gro("struttura.gro")
```

#### Da zero (a mano per test)

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

## 4. Formato config.json — Tutti i Parametri

Il file **config.json** contiene TUTTI i parametri della simulazione.
Lo schema è **stretto**: ogni chiave è obbligatoria. Per generare un template:
```python
import dmd
print(dmd.generate_config_template())
```

### 5.1 Sezione `run` — Parametri di esecuzione

```json
"run": {
    "dt": 0.001,
    "n_steps": 5000,
    "init_temperature": 300.0,
    "gen_vel": true,
    "seed": 42
}
```

| Chiave | Tipo | Obbligatorio | Descrizione |
|---|---|---|---|
| `dt` | float | sì | Passo di integrazione temporale in **picosecondi**. Valori tipici: 0.002 ps (2 fs) per sistemi con vincoli CH; 0.001 ps per sistemi tutti-atomo senza vincoli. |
| `n_steps` | int | sì | Numero totale di passi di integrazione. Il tempo totale simulato = dt × n_steps. |
| `init_temperature` | float | sì | Temperatura di generazione velocità iniziali in **Kelvin**. Usata solo se `gen_vel = true`. |
| `gen_vel` | bool | sì | Se `true`, genera velocità iniziali da una distribuzione Maxwell-Boltzmann alla temperatura `init_temperature`. Se `false`, le velocità partono da zero (utile per minimizzazione o restart). |
| `seed` | int | sì | Seme del generatore di numeri casuali. Stesso seed = stessa sequenza di velocità iniziali. Cambiando seed si ottengono diverse condizioni iniziali. |

### 5.2 Sezione `output` — Output e Traiettoria

```json
"output": {
    "trajectory_path": "traiettoria.h5",
    "nstxout": 100,
    "nstvout": 0,
    "nstenergy": 10,
    "energy_path": "energie.h5",
    "checkpoint_interval": 500,
    "checkpoint_path": "checkpoint.bin"
}
```

| Chiave | Tipo | Obbligatorio | Descrizione |
|---|---|---|---|
| `trajectory_path` | string | sì | Percorso del file H5MD per la traiettoria. Stringa vuota `""` = nessun output traiettoria. |
| `nstxout` | int | sì | Intervallo di passi per scrivere le **posizioni** nel file H5MD. 0 = non scrivere posizioni. Tipico: 50–1000. |
| `nstvout` | int | sì | Intervallo di passi per scrivere le **velocità** nel file H5MD. 0 = non scrivere velocità (risparmia spazio su disco). |
| `nstenergy` | int | sì | Intervallo di passi per scrivere le **energie** in un file H5MD separato. 0 = non scrivere. |
| `energy_path` | string | sì | Percorso del file H5MD per le energie. Usato se `nstenergy > 0`. |
| `checkpoint_interval` | int | sì | Intervallo di passi tra checkpoint (salvataggio stato). 0 = disabilitato. |
| `checkpoint_path` | string | sì | Percorso del file di checkpoint binario. Usato se `checkpoint_interval > 0` e per il checkpoint finale. |

**Relazione tra nstxout e nstvout**: le posizioni e le velocità sono
scritte in **stessi frame** (stesso passo). Se `nstxout = 50` e
`nstvout = 50`, ogni frame contiene entrambe. Se `nstvout = 0`, le
velocità non vengono scritte. Non è possibile avere frequenze diverse
per posizioni e velocità — l'intervallo effettivo è determinato da
`nstxout`; se `nstvout > 0` e diverso, il minimo comune multiplo
determina i frame effettivi (nella pratica teneteli uguali).

### 5.3 Sezione `lj` — Lennard-Jones

```json
"lj": {
    "cutoff": 1.2
}
```

| Chiave | Tipo | Obbligatorio | Descrizione |
|---|---|---|---|
| `cutoff` | float | sì | Raggio di troncamento per le interazioni LJ in **nm**. Tipico: 1.0–1.4 nm. Oltre questa distanza il contributo LJ è considerato zero. |

La lista di Verlet (lista dei vicini entro cutoff + skin) viene ricostruita
ogni 10 passi di integrazione.

### 5.4 Sezione `electrostatics` — Elettrostatica

```json
"electrostatics": {
    "coulomb_type": "pme",
    "cutoff": 1.2,
    "pme_order": 4,
    "pme_grid_spacing": 0.12,
    "ewald_coeff": 3.5
}
```

| Chiave | Tipo | Obbligatorio | Descrizione |
|---|---|---|---|
| `coulomb_type` | string | sì | Tipo di trattamento elettrostatico. Valori: `"cutoff"` (troncamento semplice), `"direct"` (Ewald direct sum con erfc), `"pme"` (Particle Mesh Ewald completo). |
| `cutoff` | float | sì | Raggio di troncamento per la parte **real-space** di Ewald in nm. Deve essere uguale o maggiore del cutoff LJ. |
| `pme_order` | int | sì | Ordine delle B-spline per l'interpolazione su griglia PME. Tipico: 4 (cubico) o 5 (più accurato, più costoso). |
| `pme_grid_spacing` | float | sì | Spaziatura desiderata della griglia PME in nm. Il numero effettivo di punti griglia viene calcolato come `ceil(L / spacing)` dove L è la dimensione della scatola. Tipico: 0.10–0.15 nm. |
| `ewald_coeff` | float | sì | Coefficiente di split di Ewald α in **nm⁻¹**. Valori tipici: 3.0–4.0. α grande → più lavoro in real-space (converge più veloce), meno in reciproco. α piccolo → più lavoro in FFT. La scelta ottimale bilancia i due costi. |

**Scelta del `coulomb_type`:**
- `"cutoff"`: troncamento netto a `cutoff`. **Attenzione**: produce artefatti
  per sistemi carichi. Usare solo per test o sistemi neutri senza cariche.
- `"direct"`: somma diretta di Ewald con termine erfc/r. Include la correzione
  per coppie escluse (1-2, 1-3). Più accurato di cutoff ma complessità O(N²).
- `"pme"`: Particle Mesh Ewald. Separa il potenziale in parte reale (erfc/r
  entro cutoff) + parte reciproca (FFT 3D su griglia) + self-energy.
  Complessità O(N log N). **Raccomandato** per sistemi con > 2000 atomi.

**Il coeffiente di Ewald α** determina quanto della forza viene calcolato
in spazio reale vs reciproco. Un valore ragionevole è `α = 3.5 / cutoff`.
A parità di cutoff, α più alto richiede più coppie in real-space ma una
griglia FFT più rada.

### 5.5 Sezione `thermostat` — Termostato

```json
"thermostat": {
    "type": "nose-hoover",
    "temperature": 300.0,
    "tau": 0.1,
    "frequency": 10.0
}
```

| Chiave | Tipo | Obbligatorio | Descrizione |
|---|---|---|---|
| `type` | string | sì | Tipo di termostato. Valori: `"none"`, `"berendsen"`, `"nose-hoover"`, `"andersen"`. |
| `temperature` | float | sì | Temperatura target in **Kelvin**. |
| `tau` | float | sì | Costante di accoppiamento in **ps** per Berendsen e Nose-Hoover. Più piccolo = accoppiamento più forte. |
| `frequency` | float | sì | Frequenza di collisione in **ps⁻¹** per Andersen. Numero medio di collisioni per atomo al picosecondo. |

**Dettaglio termostati:**

| Tipo | Meccanismo | Ensemble esatto? | Uso tipico |
|---|---|---|---|
| `"none"` | Nessun termostato (NVE) | Sì (microcanonico) | Test di conservazione energia, sistemi isolati |
| `"berendsen"` | Rescaling proporzionale delle velocità | No | Equilibrazione rapida (τ piccolo = 0.1 ps), **mai in produzione** |
| `"nose-hoover"` | Variabile estesa ζ (frizione) | Sì (canonico NVT) | Produzione (τ = 0.5–2.0 ps) |
| `"andersen"` | Collisioni stocastiche con bagno termico | Sì (canonico NVT) | Sistemi piccoli, quando si accetta dinamica non continua |

**Nose-Hoover**: è il termostato canonico esatto raccomandato per la produzione.
Introduce un grado di libertà aggiuntivo ζ (coefficiente di frizione) con massa
`Q = N_f · k_B · T₀ · τ²`. La dinamica è deterministica e reversibile.

**Berendsen**: scala le velocità come `λ = √(1 + Δt/τ · (T₀/T − 1))`.
Non produce un vero ensemble canonico (sopprime le fluttuazioni di temperatura).
Usare solo per equilibrazione iniziale.

**Andersen**: a ogni passo, ogni atomo ha probabilità `P = ν · Δt` di subire
una collisione che rimpiazza la sua velocità con un campione Maxwell-Boltzmann.
Produce un vero ensemble canonico ma distrugge la dinamica (le velocità non
sono continue).

### 5.6 Sezione `barostat` — Barostato

```json
"barostat": {
    "type": "none",
    "pressure": 1.0,
    "tau": 1.0,
    "compressibility": 4.5e-5
}
```

| Chiave | Tipo | Obbligatorio | Descrizione |
|---|---|---|---|
| `type` | string | sì | Tipo di barostato. Valori: `"none"`, `"berendsen"`, `"andersen"`. |
| `pressure` | float | sì | Pressione target in **bar**. |
| `tau` | float | sì | Costante di accoppiamento in **ps** per il barostato. |
| `compressibility` | float | sì | Comprimibilità isoterma in **bar⁻¹** (solo Berendsen). Per acqua a 300K: 4.5×10⁻⁵ bar⁻¹. |

**Dettaglio barostati:**

| Tipo | Meccanismo | Uso tipico |
|---|---|---|
| `"none"` | Volume fisso (NVE o NVT) | Simulazioni a volume costante |
| `"berendsen"` | Scaling proporzionale del volume | Equilibrazione, **non esatto** |
| `"andersen"` | Pistone esteso (grado libertà aggiuntivo ϵ) | Produzione NPT |

**Berendsen barostat**: scala il volume come `μ = 1 − Δt/τ_P · β · (P₀ − P)`.
Fattore di scala per posizioni e scatola: `μ^(1/3)`. Non riproduce fluttuazioni
di volume corrette nell'ensemble NPT.

**Andersen barostat**: usa una lagrangiana estesa con variabile ϵ (pistone)
di massa `W = N_f · k_B · T₀ · τ²`. Le posizioni scalano come `exp(ϵ̇ · Δt)`.
Più accurato del Berendsen per ensemble NPT.

### 5.7 Sezione `constraints` — Vincoli

```json
"constraints": {
    "type": "none",
    "tolerance": 1e-6
}
```

| Chiave | Tipo | Obbligatorio | Descrizione |
|---|---|---|---|
| `type` | string | sì | Tipo di vincolo. Attualmente solo `"none"`. |
| `tolerance` | float | sì | Tolleranza per il costruttore di vincoli (usata da SHAKE). |

Al momento i vincoli SHAKE sono implementati in C++ ma non ancora collegati
al SimulationEngine. I legami con idrogeno vanno attualmente integrati con
passo piccolo (dt = 0.001 ps) o congelati.

### 5.8 Esempi completi

#### NVE (microcanonico)

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

#### NVT con Nose-Hoover (produzione)

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

#### NPT con Berendsen (equilibrazione)

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

## 5. Eseguire una Simulazione — End-to-End

### 6.1 Flusso completo

```
┌─────────────────┐
│   system.json   │  ← da PDB+PSF+FF, o da GRO, o scritto a mano
└────────┬────────┘
         ↓
┌─────────────────┐
│  dmd.run()      │  ← legge system.json + config.json
│  SystemBuilder  │     costruisce SimulationConfig
└────────┬────────┘
         ↓
┌─────────────────┐
│  Engine.run()   │  ← ciclo MD: forze → integrazione → termostato → output
│                  │     scrive H5MD + checkpoint
└────────┬────────┘
         ↓
┌─────────────────┐
│  Risultati      │  → engine.positions, .velocities, .potential_energy
│  H5MD file      │  → traiettoria per analisi (MDAnalysis, VMD, Python)
│  Checkpoint     │  → restart
└─────────────────┘
```

### 6.2 Script minimale di esempio

```python
import dmd

# Carica sistema
with open("system.json") as f:
    import json
    system = json.load(f)

# Carica configurazione
with open("config.json") as f:
    config = json.load(f)

# Esegui simulazione
engine = dmd.run(system_data=system, config_data=config)

# Risultati
print(f"Passi: {engine.current_step}")
print(f"Tempo: {engine.current_time:.4f} ps")
print(f"Energia potenziale: {engine.potential_energy:.4f} kJ/mol")
```

### 6.3 Output con progress bar

```python
engine = dmd.run(
    system_data=system,
    config_data=config,
    progress_interval=100   # stampa progresso ogni 100 passi
)
```

Formato output:
```
        Step  Time (ps)    PE (kJ/mol)   Box (nm)   ns/day
--------------------------------------------------------
        500     1.0000       -85.2341     3.0000 175575.5
       1000     2.0000       -85.1234     3.0000 180234.1
```

### 6.4 Script completo PDB+PSF+FF

```python
import dmd, json

# 1. Prepara sistema
system = dmd.from_pdb("struttura.pdb")

psf = dmd.from_psf("struttura.psf", scale_14=1.0)
system["atoms"]["type"] = psf["atom_types"]
system["force_field"].update(psf["force_field"])

ff_params = dmd.load_forcefield("parametri.prm", format="charmm")
merged = dmd.merge_ff(system, ff_params)
system["force_field"].update(merged["force_field"])
system["exclusions"] = psf["exclusions"]

# 2. Configura simulazione
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

# 3. Esegui
engine = dmd.run(system_data=system, config_data=config, progress_interval=100)

# 4. Analisi base
import numpy as np
pos = np.array(engine.positions)
vel = np.array(engine.velocities)
masses = np.array(system["atoms"]["mass"])
ekin = 0.5 * np.sum(masses * np.sum(vel**2, axis=1))
kb = 0.008314462618
T = 2.0 * ekin / (3.0 * len(masses) * kb)
print(f"Temperatura finale: {T:.2f} K")
```

### 5.5 Script minimale minimizzazione + produzione

```python
import dmd

# Minimizza
info = dmd.minimize(system, n_steps=500, output={"checkpoint_path": "min.json"})

# Equilibra NVT
config_nvt = { ... thermostat: "nose-hoover" ... }
eng = dmd.run(checkpoint="min.json", config_data=config_nvt)

# Produzione NPT
config_npt = { ... thermostat: "nose-hoover", barostat: "andersen" ... }
eng2 = dmd.run(checkpoint="checkpoint.bin", config_data=config_npt)
```

### 5.6 Minimizzazione energetica

Prima di una simulazione di produzione, è buona pratica **minimizzare l'energia**
del sistema per rimuovere contatti sterici sfavorevoli:

```python
import dmd

# Minimizza
info = dmd.minimize(
    system_data=system,
    n_steps=500,           # massimo passi
    step_size=0.001,       # passo iniziale in nm
    energy_tol=0.01,       # convergenza: variazione energia < 0.01 kJ/mol per passo
    output={"checkpoint_path": "minimized_checkpoint.json"}
)

# Avvia simulazione dal minimo
engine = dmd.run(
    checkpoint="minimized_checkpoint.json",
    config_data=config
)
```

L'algoritmo è **steepest descent** con backtracking: se l'energia aumenta,
il passo viene ridotto (×0.5); se diminuisce, il passo viene aumentato (×1.2).
L'ottimizzazione si ferma quando la variazione di energia scende sotto
`energy_tol` o quando si raggiunge `n_steps`.

### 5.7 Esempio: Acqua TIP3P

L'esempio completo per una scatola d'acqua TIP3P (522 molecole, 1566 atomi)
si trova in `test_acqua/`:

```bash
cd test_acqua

# Genera sistema (PDB + PSF + PRM)
python3 generate_water_box.py

# Esegui simulazione
python3 run_water_test.py

# Visualizza traiettoria
python3 visualize_vdw.py
```

Parametri TIP3P usati:
- **O–O LJ**: σ = 0.3150 nm, ε = 0.636 kJ/mol
- **O–H bond**: k = 502416 kJ/(mol·nm²), r₀ = 0.09572 nm
- **H–O–H angle**: k_θ = 628.02 kJ/(mol·rad²), θ₀ = 104.52°
- **Cariche**: O = −0.834 e⁻, H = +0.417 e⁻
- **Box**: 2.5 nm cubica, 522 H₂O

---

## 6. Algoritmi in Dettaglio

### 6.1 Integratore Velocity Verlet

DMD usa l'integratore **Velocity Verlet**, un algoritmo simplettico (conserva
l'energia) e reversibile nel tempo:

```
v(t + ½Δt) = v(t) + (Δt/2) · a(t)              # half-kick
r(t + Δt)  = r(t) + Δt · v(t + ½Δt)            # advance
a(t + Δt)  = F(t + Δt) / m                     # forze
v(t + Δt)  = v(t + ½Δt) + (Δt/2) · a(t + Δt)  # half-kick
```

Proprietà:
- **Simplettico**: conserva l'energia totale nel limite di passo infinitesimo.
  Su passi finiti, l'energia oscilla attorno a un valore costante senza deriva.
- **Reversibile**: invertendo il segno di Δt si torna indietro.
- **Errore locale**: O(Δt³) su posizioni, O(Δt²) su velocità.
- **Errore globale**: O(Δt²) su posizioni.

Il ciclo MD completo in DMD:

```
for step in range(n_steps):
    integrator.half_kick(dt/2)              # v ← v + dt/2 · F/m
    integrator.advance(dt, box)             # r ← r + dt · v (con PBC)
    force_engine.compute(pos → F, E)        # calcola tutte le forze
    integrator.half_kick(dt/2)              # v ← v + dt/2 · F/m
    thermostat.apply(vel, dt)               # regola temperatura
    barostat.apply(pos, box, dt)            # regola pressione
    constraints.apply(pos, vel)             # vincoli (SHAKE)
    trajectory.write(pos, vel, box, step)   # H5MD
    checkpoint.save(state)                  # checkpoint periodico
```

### 6.2 Termostato Nose-Hoover

Il termostato di **Nose-Hoover** estende il sistema con un grado di libertà
aggiuntivo ζ (coefficiente di frizione). Le equazioni del moto sono:

```
ṙ_i = v_i
v̇_i = F_i/m_i − ζ · v_i
ζ̇ = (T(t) − T₀) · k_B / Q
```

dove `Q = N_f · k_B · T₀ · τ²` è la massa del termostato.

In DMD, il termostato viene applicato alla seconda metà del Velocity Verlet
con un passo di integrazione split-operator (metà prima di half-kick, metà dopo).

### 6.3 PME (Particle Mesh Ewald)

L'elettrostatica è il termine computazionalmente più costoso. La somma di
Ewald separa il potenziale Coulombiano 1/r in due parti che convergono
rapidamente:

```
1/r = erfc(αr)/r + erf(αr)/r

E_Coul = ½ Σ' q_i q_j · erfc(α r_ij) / r_ij    [real space, cutoff]
       + ½ Σ' q_i q_j · erf(α r_ij) / r_ij     [reciprocal space, FFT]
       - (k_e · α / √π) · Σ q_i²                [self-energy, sottratto]
```

**PME** (Particle Mesh Ewald) implementa la parte reciproca in tre passi:

1. **Charge spreading**: le cariche puntiformi vengono interpolate su una
   griglia 3D regolare usando B-spline di ordine n (tipicamente 4, cubiche).
   Ogni carica contribuisce a una piccola regione della griglia (n³ punti).

2. **FFT 3D**: la griglia carica viene trasformata in spazio di Fourier via
   FFT 3D (r2c half-complex), moltiplicata per la funzione di Green
   `G(k) = 4π · exp(−k²/4α²) / k² · B(k)` (dove B(k) è la correzione
   B-spline), e ritrasformata (c2r).

3. **Back-interpolation**: il potenziale reciproco su griglia viene
   re-interpolato alle posizioni atomiche usando le derivate delle B-spline,
   producendo le forze.

**Complessità**: O(N log N) via FFT, contro O(N²) della somma diretta.
Per N > 2000, PME è più veloce del direct sum a parità di accuratezza.

**Correzione esclusioni**: PME include automaticamente tutte le coppie,
anche quelle 1-2 e 1-3 che devono essere escluse. Il componente
`CoulombExclusion` sottrae il termine `q_i·q_j/r` per le coppie escluse
dopo il calcolo PME.

### 6.4 Lista di Verlet

Per le interazioni a corto raggio (LJ e Coulomb direct), DMD mantiene una
**lista di Verlet** (lista di vicini entro cutoff + skin distance).
La lista viene ricostruita ogni 10 passi per bilanciare correttezza
(gli atomi non possono percorrere più della skin distance in 10 passi)
e costo computazionale (O(N²) per ricostruire, O(N·M) per calcolare).

### 6.5 Minimizzazione Steepest Descent

```python
info = dmd.minimize(system_data, n_steps=500, step_size=0.001, energy_tol=0.01)
```

L'algoritmo:
1. Calcola forze `F = −∇E(pos)`
2. Sposta atomi lungo la forza: `x ← x + α · F/|F|`
3. Se l'energia è diminuita (`E_nuova < E_precedente`): accetta, raddoppia α
4. Se l'energia è aumentata: rifiuta, dimezza α, riprova
5. Ripeti fino a convergenza (`|ΔE| < tol`) o `n_steps`

---

## 7. Il Formato H5MD

H5MD (**HDF5 Molecular Dynamics**) è uno standard aperto per la memorizzazione
di traiettorie MD basato su HDF5. È stato sviluppato dalla comunità come
alternativa aperta ai formati proprietari o legacy.

### 8.1 Perché H5MD e non DCD, XTC, o TRR?

| Caratteristica | H5MD | DCD (CHARMM) | XTC (GROMACS) | TRR (GROMACS) |
|---|---|---|---|---|
| **Aperto** | Sì, standard documentato | No, binario proprietario | No, specifico GROMACS | No, specifico GROMACS |
| **Auto-descrittivo** | Sì (HDF5 contiene i metadati) | No | No | No |
| **Accesso parziale** | Sì (leggi solo posizioni, solo certi frame) | Sì (ma complesso) | No (decomprime) | Sì |
| **Tipi multipli** | Posizioni, velocità, forze, box, parametri | Solo posizioni | Solo posizioni | Posizioni + velocità + forze |
| **Estendibile** | Sì (dataset chunked) | Sì (append) | No (riscrive) | Sì (append) |
| **Interoperabile** | MDAnalysis, VMD, h5py, Julia, MATLAB | CHARMM, MDAnalysis | GROMACS, MDAnalysis, VMD | GROMACS, MDAnalysis |
| **Compressione** | Sì (filtri HDF5: gzip, szip, lzf) | No | Sì (proprietaria) | No |

**Vantaggi decisivi di H5MD:**
- **Portabile**: un file H5MD può essere letto da Python (h5py), VMD,
  MDAnalysis, Julia, MATLAB, e da qualsiasi linguaggio con bindings HDF5.
- **Auto-descrittivo**: i nomi dei dataset e i loro attributi (unità,
  step temporali) sono dentro il file stesso. Non serve un file di
  parametri separato.
- **Accesso casuale**: puoi leggere solo gli ultimi 100 frame di un
  file da 1 milione di frame senza decomprimere tutto.
- **Unico file**: traiettoria, box, tempo e passo sono nello stesso file.
  Le energie possono stare in un file separato o nello stesso.

### 8.2 Struttura di un file H5MD

```
traiettoria.h5
├── /h5md/
│   └── version = [1, 1]
├── /particles/
│   └── /atoms/
│       ├── /position/
│       │   ├── value     [n_frames, n_atoms, 3]  double
│       │   ├── time      [n_frames]               double
│       │   ├── step      [n_frames]               int
│       │   └── specification (attributi: units = "nm", ...)
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

I dataset `value` sono **chunked** e **extendibili**: nuovi frame vengono
aggiunti con `H5Dset_extent` + hyperslab, senza dover ricreare il file.

### 8.3 Leggere H5MD con Python

```python
import h5py
import numpy as np

with h5py.File("traj.h5", "r") as f:
    pos = f["/particles/atoms/position/value"][:]    # (n_frames, n_atoms, 3)
    step = f["/particles/atoms/position/step"][:]    # (n_frames,)
    time = f["/particles/atoms/position/time"][:]    # (n_frames,)
    box = f["/particles/atoms/box/edges/value"][:]   # (n_frames, 3, 3)

print(f"Frame: {pos.shape[0]}, Atomi: {pos.shape[1]}")
```

Con MDAnalysis:
```python
import MDAnalysis as mda
u = mda.Universe("traj.h5")
for ts in u.trajectory:
    print(ts.frame, ts.time, ts.positions.shape)
```

---

## 8. Checkpoint e Restart

### 9.1 Checkpoint in Python (JSON)

```python
import dmd
import json

# Dopo la simulazione
dmd.write_checkpoint("checkpoint.json", engine, system_data)

# Per riprendere
engine = dmd.run(checkpoint="checkpoint.json", config_data=config)
```

Il checkpoint JSON contiene:
- `format: "dmd_checkpoint"`, `version: 1`
- `system`: l'intero `system.json` (topologia + parametri), così il restart
  **non dipende dal force field originale** (il checkpoint è autocontenuto)
- `state`: passo corrente, tempo, posizioni, velocità, box_size, energia

### 9.2 Checkpoint in C++ (binario)

Durante la simulazione, se `checkpoint_interval > 0`, il motore C++ scrive
un checkpoint binario (`magic = 0x444D4450`) con:
- Stato completo: step, time, PE, posizioni, velocità, forze, masse, cariche,
  tipi atomici
- Stato della lista di Verlet: `step_since_rebuild` (garantisce restart
  **bit-exact** — la traiettoria dopo restart è identica a quella continua)

### 9.3 Workflow completo con checkpoint

```python
import dmd, json

# Passo 1: avvio nuovo
config = { ... }  # config.json dict
system = { ... }  # system.json dict
engine = dmd.run(system_data=system, config_data=config)  # produce checkpoint

# Passo 2: continua da checkpoint
engine2 = dmd.run(checkpoint="checkpoint.bin", config_data=config)

# Passo 3: seconda continuazione
engine3 = dmd.run(checkpoint="checkpoint.bin", config_data=config)
```

### 9.4 Minimizzazione + checkpoin + produzione

```python
# 1. Minimizza
info = dmd.minimize(system, n_steps=500, output={"checkpoint_path": "min.json"})

# 2. Equilibra NVT
config_nvt = { ... thermostat: nose-hoover ... }
eng = dmd.run(checkpoint="min.json", config_data=config_nvt)

# 3. Produzione NPT
config_npt = { ... thermostat: nose-hoover, barostat: andersepy ... }
eng2 = dmd.run(checkpoint="checkpoint.bin", config_data=config_npt)
```

---

## 9. Appendice: Best Practices

### Scelta del passo temporale (dt)

| Sistema | dt raccomandato | Note |
|---|---|---|
| Tutti atomi (no vincoli H) | 0.001 ps (1 fs) | Sicuro per ogni sistema |
| Con vincoli C–H, O–H, N–H | 0.002 ps (2 fs) | Standard MD biomolecolare |
| Coarse-grained | 0.005–0.020 ps | Dipende dalla risoluzione CG |

### Scelta dei cutoff

| Interazione | Cutoff raccomandato | Note |
|---|---|---|
| Lennard-Jones | 1.0–1.4 nm | 1.2 nm è lo standard biomolecolare |
| Coulomb (direct) | 1.0–1.4 nm | Deve essere uguale per LJ e Coulomb |
| PME grid spacing | 0.10–0.15 nm | Più fine = più accurato, più lento |
| PME spline order | 4 | Ordine 4 (cubico) è il miglior rapporto qualità/costo |

### Costanti di accoppiamento

| Scopo | τ termostato | τ barostato |
|---|---|---|
| Equilibrazione rapida | 0.1 ps | 1.0 ps |
| Produzione NVT (Nose-Hoover) | 1.0–2.0 ps | — |
| Produzione NPT | 1.0–2.0 ps | 5.0–10.0 ps |

### Performance

- **LJ + Coulomb direct**: O(N·M), dominano per sistemi piccoli (< 2000 atomi)
- **PME**: O(N log N), più veloce di direct sum per N > 2000
- **Verlet list**: ricostruzione ogni 10 passi, costo O(N²) ma ammortizzato
- **ns/day**: per 1566 atomi (TIP3P, 2.5 nm box, NVT, PME): ~1.2 ns/giorno
  su un core Apple Silicon M3

### Debug e validazione

1. **Conservazione energia**: in NVE, l'energia totale (E_pot + E_kin) deve
   oscillare senza deriva. Drift > 0.5% su 5000 step indica un problema.
2. **Restart bit-exact**: dopo un checkpoint, la traiettoria deve essere
   identica a quella continua. DMD garantisce questo per il formato binario.
3. **Forze non NaN**: controlla che le forze non diventino 10¹⁵ — tipico
   segno di coordinate errate (unità sbagliate, colonne scambiate).
4. **Temperatura**: dopo equilibrazione NVT, la temperatura deve oscillare
   attorno al target con fluttuazioni `√(2/N_f) · T₀` (~5% per 1000 atomi).

---

## 10. Appendice: Riepilogo Conversioni Unità

### CHARMM → DMD (eseguito automaticamente da charmm.py)

| Grandezza | CHARMM (input) | DMD (interno) | Fattore |
|---|---|---|---|
| σ (LJ) | Rmin/2 in Å | σ in nm | ÷10 × 2 ÷ 2^(1/6) |
| ε (LJ) | kcal/mol | kJ/mol | × 4.184 |
| k (bond) | kcal/(mol·Å²) | kJ/(mol·nm²) | × 4.184 × 100 |
| r₀ (bond) | Å | nm | ÷ 10 |
| k_θ (angle) | kcal/(mol·rad²) | kJ/(mol·rad²) | × 4.184 |
| θ₀ (angle) | gradi | radianti | × π/180 |
| k_φ (dihedral) | kcal/mol | kJ/mol | × 4.184 |
| φ₀ (dihedral) | gradi | radianti | × π/180 |

### PDB → DMD

Le coordinate PDB sono in **Ångstrom** (tradizionalmente). DMD le vuole
in **nm**. La conversione (÷10) è a carico dell'utente o del runner
script. Il parser `from_pdb()` lascia i valori raw.

### DMD → SI

| DMD | SI |
|---|---|
| 1 nm | 10⁻⁹ m |
| 1 kJ/mol | 1.66054 × 10⁻²¹ J per molecola |
| 1 ps | 10⁻¹² s |
| 1 kJ/(mol·nm) | 1.66054 × 10⁻¹² N per molecola |
