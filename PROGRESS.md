# Progresso sviluppo DMD

## Stato attuale
- **Fase corrente:** Fase 1/2 — Core C++ (CPU)
- **Ultimo gate superato:** **Gate Produzione (end-to-end)**
- **Ultima feature completata:** CLI + config.json strict + build_simulation orchestrator

## Milestones
- [x] **Gate Argon NVE** — simulazione Argon 256 atomi, 5000 step, energia conservata (drift < 0.5%)
- [x] **Gate checkpoint/restart** — restart produce traiettorie identiche (bit-exact)
- [x] **Gate termostati/barostati** — 3 termostati + 2 barostati, testati singolarmente
- [x] **Gate NPT cycle** — ciclo NPT (Berendsen T+P) su Argon 108 atomi, T=120K, P=100 bar
- [x] **dmdin binary I/O** — formato binario input + JSON config reader
- [x] **H5MD trajectory output** — scrittura traiettorie in formato H5MD
- [x] **Produzione** — run completo da file .dmdin + config.json a traiettoria H5MD
- [ ] **Fase 1: Core C++ completata**

---

## Cronologia

### 2026-06-16 — Avvio sviluppo

- File: `PROGRESS.md` (questo file)
- Piano di sviluppo definito e condiviso
- Skills create: `dmd-code-style`, `dmd-workflow`
- Design document: `design/force_engine.md`
- Diagramma architetturale: `design/architettura.mmd`

### 2026-06-16 — Step 1: Toolchain + Cell + SystemData + HarmonicBond + Integrator

- **Toolchain:** CMakeLists.txt radice, `src/`, `tests/` con moduli core/force/integrate. GoogleTest via FetchContent. C++20, Clang.
- **Cell:** `core/cell.h/.cpp` — 4 tipi (orthorhombic/triclinic/hexagonal/rhombic_dodecahedron), `minimum_image()`, `volume()`, `wrap()`, `scale()`
- **SystemData:** `core/system_data.h/.cpp` — allocazione SoA centralizzata, accesso via `std::span`
- **Error handling:** `core/error.h/.cpp` — `ForceException` con enum `ForceError`
- **ForceEngine:** `force/force_engine.h/.cpp` — orchestrator con `vector<unique_ptr<ForceComponent>>`
- **ForceComponent:** `force/force_component.h/.cpp` — interfaccia astratta
- **HarmonicBond:** `force/harmonic_bond.h/.cpp` — E = ½k(r−r₀)², forze con minimum_image
- **HarmonicAngle:** `force/harmonic_angle.h/.cpp` — E = ½k_θ(θ−θ₀)²
- **PeriodicDihedral:** `force/periodic_dihedral.h/.cpp` — E = k_φ[1+cos(nφ−φ₀)], forze implementate con gradienti analitici
- **LennardJones:** `force/lennard_jones.h/.cpp` — LJ E = 4ε[(σ/r)¹²−(σ/r)⁶], Verlet list
- **CoulombDirect:** `force/coulomb_direct.h/.cpp` — Ewald direct sum, erfc/r
- **CoulombPME:** `force/coulomb_pme.h/.cpp` — stub (FFTW pending)
- **Integrator:** `integrate/integrator.h/.cpp` — Velocity Verlet split in half_kick + advance
- **Constraints:** `integrate/constraints.h/.cpp` — SHAKE iterativo
- **Test:** 11 test passati (10 unit + 1 integration)

### 2026-06-16 — Gate Argon NVE + fix bug LJ

- **Bug LJ (doppio):** formula divideva per `r` invece di `r²`, segno forze invertito su i/j
- **Effetto:** atomi attratti invece di respinti → collasso → energia esplosa a 10⁴⁵
- **Fix:** `f_scalar = 24ε(2sr⁶²−sr⁶)/r²`, accumulation `i -= / j +=`
- **Gate:** test integratione Argon NVE con 256 atomi FCC, densità 1.4 g/cm³, T=85K
- **Drift:** 0.15% su 5000 step, tolleranza 0.5%
- **Stato:** 11/11 test passano

### 2026-06-16 — CoulombPME con FFTW (completamento Fase 2)

- **FFTW:** integrato via `PkgConfig::FFTW` (Homebrew 3.3.11 su Apple Silicon)
- **CoulombPME:** implementazione completa SPME (Particle Mesh Ewald)
  - **B-spline:** cubic cardinal B-spline (ordine 4), supporto 4 punti per dimensione
  - **Charge spreading:** interpolazione B-spline su griglia 3D con PBC
  - **FFT (r2c/c2r):** FFTW 3D in-place per efficienza
  - **Green's function:** 4πk_e/V * exp(−k²/4α²)/k² * B(k) (con correzione B-spline)
  - **Peso Hermitiano:** w=2 per 0<kx<nx/2 (r2c dimezzato)
  - **Self-energy:** −k_e·α/√π·Σq²
  - **Forze:** interpolazione dal potenziale reciproco con derivate B-spline e fattori dimensioni nx/Lx, ny/Ly, nz/Lz
- **Test:** 6 test (B-spline, self-energy, simmetria forze, conservazione energia NVE breve)
- **Stato:** 12/12 test passano (10 unit + 1 integration + 1 PME)

### 2026-06-16 — Simulation Engine + checkpoint/restart (Gate checkpoint)

- **SimulationEngine:** classe orchestrazione con `step()`, `run()`, `save_checkpoint()`, `load_checkpoint()`
- **SimulationConfig:** struct dati + `build_simulation()` factory function
- **Checkpoint:** formato binario (magic 0x444D4450, version 1), stato completo system + LJ Verlet
- **ForceEngine:** esposto `components()` per serializzazione checkpoint
- **LennardJones:** esposto `step_since_rebuild()` / `set_step_since_rebuild()` per restart deterministico
- **Bug forces (doppio):** forze non azzerate all'inizio di `step()` causavano accumulo sul passo precedente → forze raddoppiate sulla prima chiamata `compute()`.
- **Bug off-by-one:** `sys.step` impostato prima dell'incremento di `step()` → checkpoint salvava step-1.
- **Gate checkpoint/restart:** restart da checkpoint produce traiettoria bit-exact (32 atomi Ar, 200 step, split 100+100)
- **Nuovi file:** `src/sim/` (engine, config, checkpoint, CMakeLists), `tests/unit/sim/`, `tests/integration/test_checkpoint_restart.cpp`
- **Stato:** 14/14 test passano (11 unit + 3 integration)

### 2026-06-16 — Thermostat + Barostat modules (Gate termostati/barostati)

- **Thermostat base:** `thermostat/thermostat.h/.cpp` — classe astratta + funzioni `kinetic_energy()`, `temperature()`
- **AndersenThermostat:** collisioni stocastiche con bagno termico a frequenza ν
- **BerendsenThermostat:** rescaling velocities con τ di accoppiamento
- **NoseHooverThermostat:** friction coefficient ζ, massa termostatica Q = N_f·k_B·T₀·τ²
- **Barostat base:** `barostat/barostat.h/.cpp` — classe astratta + funzioni `virial()`, `pressure()`
- **BerendsenBarostat:** scaling isotropico posizioni/cella, compressibilità isothermal
- **AndersenBarostat:** pistone esteso (ϵ̇, W), equazioni di Langevin per volume
- **Integrazione:** SimulationEngine accetta Thermostat/Barostat via costruttore o setter, li applica dopo la seconda metà del Velocity Verlet
- **Costanti fisiche:** `core/constants.h` con KB, PI
- **Nuove directory:** `src/thermostat/`, `src/barostat/`, `tests/unit/thermostat/`, `tests/unit/barostat/`
- **Stato:** 16/16 test passano (13 unit + 3 integration)

### 2026-06-16 — dmdin binary format + reader/writer (System Layer input)

- **Formato dmdin:** formato binario per input sistema (magic 0x444D444E, version 1)
  - Header 120B: MAGIC·VERSION·n_atoms·cell_type·n_types·pos_flags·offset_ff·offset_pos·box[9]
  - System section: per-type masses/charges
  - FF section: LJ params, bonded terms (bonds, angles, dihedrals)
  - Positions section: pos/vel (float→double), atom_types (int32)
- **Funzioni:** `read_dmdin(path)` → cfg, `write_dmdin(path, cfg)`, `apply_json_config(cfg, json)`
- **JSON config:** parametri runtime (dt, n_steps, cutoff) via nlohmann/json (FetchContent)
- **SimulationConfig:** aggiunto `atom_types` per-atomo
- **build_simulation:** ora usa `cfg.atom_types` se presente
- **Nuova dipendenza:** nlohmann/json v3.11.3 (header-only)
- **Nuovo modulo:** `src/sysbin/` + `tests/unit/sysbin/`
- **Stato:** 18/18 test passano (14 unit + 4 integration)

### 2026-06-16 — H5MD trajectory output

- **H5MDWriter:** classe in `src/trajectory/h5md_writer.h/.cpp`
  - Formato H5MD standard: `/h5md/version = [1,1]`, `/particles/atoms/{position,velocity,box/edges}`
  - Dataset estendibili (chunked, first dim = H5S_UNLIMITED) per append di frame
  - Scrive: posizioni (x,y,z interleaved), velocità (opzionale), box 3×3, time+step per frame
- **Integrazione SimulationEngine:** `Config::trajectory_path`, `Config::trajectory_interval`
  - `run()` crea H5MDWriter se `trajectory_interval > 0`
  - `save_frame()` chiama `writer.write_frame(sys, cell_);`
  - `build_simulation()` propaga `cfg.trajectory_path` → `engine_cfg.trajectory_path`
- **Nuova dipendenza:** HDF5 2.1.1 (Homebrew, via `find_package(HDF5 COMPONENTS C)`)
- **Nuovo modulo:** `src/trajectory/` + `tests/unit/trajectory/`
- **Stato:** 19/19 test passano (15 unit + 4 integration)

### 2026-06-16 — Gate Produzione (end-to-end: .dmdin + config.json → H5MD)

- **SimulationConfig esteso:**
  - `ThermostatConfig` + `BarostatConfig` struct per coupling params
  - `nstxout`, `nstvout`, `nstenergy` separati (→ sostituisce `trajectory_interval`)
  - `energy_path`, `constraint_type`, `constraint_tolerance`, `init_temperature`, `gen_vel`, `seed`
  - `coulomb_type`, `coulomb_cutoff`, `pme_order`, `pme_grid_spacing`, `ewald_coeff`
- **JSON config strict:** `src/sim/json_config.h/.cpp`
  - Validazione contro schema: tutte le chiavi richieste, chiavi sconosciute → errore
  - `generate_config_template()` per `dmd config --template`
- **build_force_engine:** ora crea CoulombDirect/PME in base a `coulomb_type`
- **build_simulation:** crea termostato/barostato da config + genera velocità Maxwell-Boltzmann se `gen_vel=true`
- **CLI:** `src/main/main.cpp` — `dmd run <system.dmdin> <config.json>` + `dmd config --template`
- **Gate Produzione:** `tests/integration/test_production.cpp` — Ar FCC, write dmdin + config.json → read → run → verifica H5MD
- **Stato:** 20/20 test passano (15 unit + 5 integration)

### 2026-06-16 — Gate NPT cycle (Berendsen T+P su Argon liquido)

- **Gate test:** `tests/integration/test_npt_berendsen.cpp`
  - 108 Ar atomi FCC, densità iniziale 1.0 g/cm³
  - Berendsen thermostat (T=120K, τ=0.1 ps) + Berendsen barostat (P=100 bar, τ=1.0 ps)
  - 2000 step, dt=0.002 ps
- **Risultati:** temperatura entro 10% del target, densità convergente verso liquido Ar
- **Stato:** 17/17 test passano (13 unit + 4 integration)

### 2026-06-16 — Python API Layer (Fase 6)

- **pybind11 module** `python/dmd/core.cpp` → `_dmd_core.cpython-314-darwin.so`
  - `SimulationConfig` binding with scalar fields + numpy array properties (mass, charges, positions, velocities, atom_types, LJ params, bond params)
  - `EngineWrapper` class with `run()`, `current_step`, `current_time`, `positions`, `velocities`, `potential_energy` (as numpy arrays)
  - Free functions: `read_dmdin`, `write_dmdin`, `apply_json_config`, `generate_config_template`, `build_cfg`
- **Python package** `python/dmd/`:
  - `__init__.py` — public API exports
  - `system.py` — `SystemBuilder` class (system.json → SimulationConfig)
  - `runner.py` — `run()` orchestrator (system.json/config.json → Engine)
  - `convert/pdb.py` — PDB parser
  - `convert/psf.py` — PSF parser (bonds, angles, dihedrals, impropers)
  - `convert/gro.py` — GRO parser
  - `convert/__init__.py` — re-exports
- `pyproject.toml` for pip install
- **CMakeLists.txt** builds `_dmd_core.so` into `python/dmd/`
- **Stato:** 20/20 test C++ passano + Python bindings verificati manualmente (build, import, config, system, run)

### 2026-06-17 — ForceField parser (design fase 7)

- **Architettura force field generica** progettata:
  - `python/dmd/forcefield/` package con `base.py`, `registry.py`, `charmm.py`
  - Classe astratta `ForceFieldParser` con metodo `parse(path) → FFParams`
  - `FFParams` struttura: atom_types, bonds, angles, dihedrals, impropers
  - `merger.py` — abbina connettività PSF + parametri FF per tipo atomico
  - Supporto futuro: AMBER (`amber.py`), OPLS (`opls.py`), GROMOS
- **Workflow utente:** `from_pdb + from_psf + load_forcefield + merge_ff → dmd.run()`
- **Da implementare:** parser CHARMM PAR/PRM, merger, integrazione in __init__.py

### 2026-06-17 — ForceField parser (implementazione fase 7)

- **File creati e modificati:**
  - `python/dmd/forcefield/__init__.py` — public API exports
  - `python/dmd/forcefield/base.py` — `FFParams` dataclass (atom_types, bonds, angles, dihedrals, impropers)
  - `python/dmd/forcefield/registry.py` — `register()` + `load_forcefield(format)` dispatcher
  - `python/dmd/forcefield/charmm.py` — CHARMM PAR/PRM parser con wildcard `X`, multi-term dihedrals
  - `python/dmd/forcefield/merger.py` — `merge_ff(psf_data, ff_params)` → abbina connettività PSF + parametri FF per tipo atomico + matrice LJ combinata (Lorentz-Berthelot)
  - `python/dmd/convert/psf.py` — aggiornato per estrarre anche `atom_types` dalla sezione !NATOM
  - `python/dmd/__init__.py` — esporta `load_forcefield`, `merge_ff`
- **Flusso di lavoro:** `from_pdb + from_psf + load_forcefield + merge_ff → system → dmd.run()`
- **Test:** end-to-end con PDB (3 atomi) + PSF + CHARMM PAR → build → run → 50 step
- **Problemi risolti:**
  - Multi-term dihedrals: due entry con (i,j,k,l) uguale ma periodo/k_phi/phi0 diverso
  - Wildcard `X` per impropers
  - Tipi atomici: PDB e PSF possono avere nomi diversi → l'utente deve sovrascrivere `atoms.type` con `psf['atom_types']`
- **Stato:** 20/20 test C++ passano + forcefield pipeline verificata manualmente

### 2026-06-17 — Minimizer + Checkpoint/Continuation (fase 7.2)

- **File creati e modificati:**
  - `python/dmd/checkpoint.py` — `write_checkpoint(path, engine, system_data)`, `read_checkpoint(path)`, `checkpoint_to_system_data(ckpt)`
  - `python/dmd/minimizer.py` — `minimize(system_data, n_steps, step_size, energy_tol, output)` steepest descent con backtracking
  - `python/dmd/runner.py` — `run()` ora accetta `checkpoint=` per continuation
  - `python/dmd/core.cpp` — esposta `forces` come `ndarray (n_atoms, 3)` su EngineWrapper
  - `python/dmd/system.py` — `use_lj=True` e `lj_cutoff` automatici quando ci sono parametri LJ
  - `python/dmd/__init__.py` — esporta `write_checkpoint`, `read_checkpoint`, `minimize`
- **Formato checkpoint.json:** autocontenuto (system + state), separato da config.json
- **Flusso:** minimize → checkpoint → run(checkpoint=, config_data=) → checkpoint → run(checkpoint=, ...)
- **Test:** fresh run (50), minimize (50) → run (50), checkpoint after MD (50) → continuation (100) → second continuation (200)
- **Stato:** 20/20 test C++ passano + minimizer/checkpoint pipeline verificata

### 2026-06-17 — Status reporter + run_n (fase 7.3)

- **File creati e modificati:**
  - `python/dmd/status.py` — `SimulationStatus` class: timer, progress table stile GROMACS (step, time, PE, box, ns/day), `start()/step()/finish()` lifecycle
  - `python/dmd/runner.py` — `run()` ora ha parametro `progress_interval=100`: chunked execution con `engine.run_n(chunk)` + status live tra i chunk
  - `python/dmd/core.cpp` — aggiunte: `run_n(steps)`, `box_size`, `n_atoms` su `EngineWrapper`
- **Formato output:**
  ```
        Step  Time (ps)    PE (kJ/mol)   Box (nm)   ns/day
  --------------------------------------------------------
        500     1.0000        -0.0031     3.0000 175575.5
        ...
    [5000 steps in 2.3s]  ns/day=178507.4
  ```
- **Test:** minimize (100) → NVE (200) → continuation NVE (500) con progress; 20/20 C++ test passano

### 2026-06-17 — Phase 0 Pre-Test Acqua: PME exclusion, FFTW threading, PSF exclusions (FIX)

- **Fase 0.3 — FFTW3 multi-thread:**
  - `src/force/coulomb_pme.cpp`: chiama `fftw_init_threads()` + `fftw_plan_with_nthreads(omp_get_max_threads())` nel costruttore
  - `src/force/CMakeLists.txt`: link `libfftw3_omp` + `libomp` (Homebrew su Apple Silicon)
- **Fase 0.2 — PSF exclusions (Python):**
  - `python/dmd/convert/psf.py`: nuova funzione `_exclusion_pairs(bonds, angles, dihedrals, scale_14)` che genera lista coppie escluse 1-2, 1-3
  - `from_psf()` ora accetta `scale_14` e restituisce `exclusions` nel dict
  - `python/dmd/system.py`: `SystemBuilder.from_system_json()` legge `exclusions` e popola `cfg.excl_i/excl_j`
- **Fase 0.1 — PME exclusion correction (C++):**
  - Nuovo componente `src/force/coulomb_exclusion.h/.cpp`: `CoulombExclusion` sottrae `q_i·q_j/r` per coppie escluse (corregge PME che include tutte le coppie)
  - `build_force_engine()`: quando `coulomb_type == "pme"`, aggiunge **entrambi** `CoulombDirect` (real-space erfc) + `CoulombPME` (reciproco FFT) + `CoulombExclusion` (correzione)
  - Quando `coulomb_type == "direct"`, aggiunge `CoulombDirect` + `CoulombExclusion` (se ci sono esclusioni)
  - `SystemData` + `SimulationConfig`: aggiunti `excl_i, excl_j` (flat pair arrays)
  - `python/dmd/core.cpp` (pybind11): esposti `excl_i`, `excl_j` su `SimulationConfig`
- **Bug critici risolti:**
  - CHARMM parser: aggiunta conversione unità LJ (kcal→kJ, A→nm, deg→rad) — sigma era in A, epsilon in kcal
  - Merger: condizionale `if sig_i and sig_j` falliva con sigma=0 (H) → sigma=1.0 per H-OH2. Fix: calcola sempre media aritmetica, condizionale solo su epsilon.
  - SystemBuilder: `pos[:,0]` crea vista numpy non contigua → pybind setter legge dati contigui ignorando stride → coordinate lette in colonne errate (X→Z shift). Fix: `.copy()` su ogni slice colonna.
  - H5MD non scritto con `engine.run_n()` (usato per progress reporting): `run_n()` chiamava `engine_.step()` bypassando `run()` che crea H5MDWriter. Fix: `EngineWrapper::run_n()` ora crea writer e salva frame a intervallo configurato.
- **Test:** 20/20 test C++ passano + water test 522 TIP3P (1566 atomi) NVT nose-hoover 300K:
  - 500 step: PE~1661 kJ/mol (1.06/atom), stabile, forces O(1000) kJ/(mol·nm)
  - 2000 step (40 frame con nstxout=50): PE stabile ~1400-2000 kJ/mol, 138s, ns/day=1.3
  - H5MD (3MB, 40 frame) e .npy generati in test_acqua/
  - HDF5 file lock issue su macOS con `H5F_ACC_TRUNC` dopo run consecutivi — risolto pulendo file pre-esistenti
- **Visualizzazione:** script `test_acqua/visualize_vdw.py` → MP4 con O (rosso) e H (bianco/grigio), Z-sorted, VDW-like scatter
- **Infrastruttura:** `.gitignore` aggiornato con `test_acqua/*`
- **Gate:** —
