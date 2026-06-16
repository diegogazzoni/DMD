# Progresso sviluppo DMD

## Stato attuale
- **Fase corrente:** Fase 1/2 — Core C++ (CPU)
- **Ultimo gate superato:** nessuno
- **Ultima feature completata:** — (nessuna, avvio sviluppo)

## Milestones
- [ ] **Gate Argon NVE** — simulazione Argon 500 atomi, 10⁶ step, energia conservata

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
- **PeriodicDihedral:** `force/periodic_dihedral.h/.cpp` — E = k_φ[1+cos(nφ−φ₀)] (energy only, forces stub)
- **LennardJones:** `force/lennard_jones.h/.cpp` — LJ E = 4ε[(σ/r)¹²−(σ/r)⁶], Verlet list
- **CoulombDirect:** `force/coulomb_direct.h/.cpp` — Ewald direct sum, erfc/r
- **CoulombPME:** `force/coulomb_pme.h/.cpp` — stub (FFTW pending)
- **Integrator:** `integrate/integrator.h/.cpp` — Velocity Verlet split in half_kick + advance
- **Constraints:** `integrate/constraints.h/.cpp` — SHAKE iterativo
- **Test:** 7 unit test passati (Cell, SystemData, HarmonicBond, LennardJones, Coulomb, ForceEngine, Integrator)
- **Gate:** — (nessuno, prerequisiti del gate Argon)
