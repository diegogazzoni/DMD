# Progresso sviluppo DMD

## Stato attuale
- **Fase corrente:** Fase 1/2 — Core C++ (CPU)
- **Ultimo gate superato:** **Gate NPT cycle**
- **Ultima feature completata:** NPT Berendsen integration test su Argon liquido

## Milestones
- [x] **Gate Argon NVE** — simulazione Argon 256 atomi, 5000 step, energia conservata (drift < 0.5%)
- [x] **Gate checkpoint/restart** — restart produce traiettorie identiche (bit-exact)
- [x] **Gate termostati/barostati** — 3 termostati + 2 barostati, testati singolarmente
- [x] **Gate NPT cycle** — ciclo NPT (Berendsen T+P) su Argon 108 atomi, T=120K, P=100 bar
- [ ] **I/O PDB/GRO** — lettura/scrittura topology e coordinate
- [ ] **Produzione** — run completo con output traiettoria ed energia

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

### 2026-06-16 — Gate NPT cycle (Berendsen T+P su Argon liquido)

- **Gate test:** `tests/integration/test_npt_berendsen.cpp`
  - 108 Ar atomi FCC, densità iniziale 1.0 g/cm³
  - Berendsen thermostat (T=120K, τ=0.1 ps) + Berendsen barostat (P=100 bar, τ=1.0 ps)
  - 2000 step, dt=0.002 ps
- **Risultati:** temperatura entro 10% del target, densità convergente verso liquido Ar
- **Stato:** 17/17 test passano (13 unit + 4 integration)
