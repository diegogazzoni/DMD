# Force Engine — Specifica tecnica

Data: 2026-06-16
Documento architetturale primario: ARCHITECTURE.md
Questo documento è parte della definizione della Fase 2 del progetto DMD.

---

## 1. Architettura a componenti

Il Force Engine scompone il calcolo delle forze in componenti indipendenti,
ciascuno responsabile di un termine energetico. Ogni componente implementa
l'interfaccia `ForceComponent`:

```cpp
class ForceComponent {
public:
    virtual ~ForceComponent() = default;
    virtual std::string name() const = 0;
    virtual void compute(
        std::span<const double> pos_x,
        std::span<const double> pos_y,
        std::span<const double> pos_z,
        const Cell& cell,
        std::span<double> forces_x,
        std::span<double> forces_y,
        std::span<double> forces_z,
        double& energy_out
    ) = 0;
};
```

### Lista fissa di componenti

Registrata una sola volta all'inizializzazione, invariante per l'intera simulazione:

| # | Componente | Categoria | Dipende da |
|---|------------|-----------|------------|
| 1 | `HarmonicBond` | Bonded | spanning tree topologia |
| 2 | `HarmonicAngle` | Bonded | angoli |
| 3 | `PeriodicDihedral` | Bonded | diedri (proper + improper) |
| 4 | `LennardJones` | Non-bonded (short) | pair list + lookup table |
| 5 | `CoulombDirect` | Non-bonded (short) | stessa pair list di LJ |
| 6 | `CoulombPME` | Long-range | griglia 3D + FFT |
| 7 | `ExternalForce` | SimulationModule | campo esterno (es. elettrico) |
| 8 | `Constraints` | Vincoli | SHAKE/RATTLE, eseguito post-integrazione |

### Energy tracking

Ogni componente riporta la propria energia in un double separato. Il ForceEngine
accumula in un array `component_energies[8]` esposto al SimulationEngine per
debug e logging.

```cpp
struct EnergyReport {
    std::array<double, 8> per_component;
    double total;
};
```

---

## 2. ForceEngine

Classe orchestrator che possiede e coordina i componenti.

```cpp
class ForceEngine {
public:
    void add_component(std::unique_ptr<ForceComponent> component);
    void set_backend(std::unique_ptr<ForceEngineBackend> backend);

    EnergyReport compute(
        SystemData& sysdata,
        const Cell& cell
    );

private:
    std::vector<std::unique_ptr<ForceComponent>> components_;
    std::unique_ptr<ForceEngineBackend> backend_;
};
```

### Ordine di esecuzione

Fisso, dichiarato all'inizializzazione:

```
 0. HarmonicBond
 1. HarmonicAngle
 2. PeriodicDihedral
 3. LennardJones (+ CoulombDirect, pair list condivisa)
 4. CoulombPME
 5. ExternalForce  (se registrato da SimulationModule)
---- integratore ----
 6. Constraints     (post-integrazione, modifica pos/vel)
```

Nota: LennardJones e CoulombDirect condividono la stessa Verlet list per
efficienza — il loop delle coppie calcola entrambi i contributi
contemporaneamente.

---

## 3. SystemData (centralizzato)

Unico contenitore di tutta la memoria dello stato di simulazione.

```cpp
struct SystemData {
    size_t n_atoms;

    // Posizioni e velocità (input/output)
    std::vector<double> pos_x, pos_y, pos_z;
    std::vector<double> vel_x, vel_y, vel_z;

    // Forze (output del ForceEngine, input dell'Integrator)
    std::vector<double> forces_x, forces_y, forces_z;

    // Proprietà atomiche (costanti)
    std::vector<double> masses;
    std::vector<double> charges;
    std::vector<int>    atom_types;

    // Energia potenziale (output)
    double potential_energy;

    // Passo corrente e tempo
    size_t step;
    double time;
};
```

### Layout SoA

I vettori sono allocati contiguamente in heap standard. I componenti
ricevono `std::span<double>` ai piani X/Y/Z, non l'intero SystemData.

### Ownership

Il `SimulationEngine` possiede il SystemData. Il ForceEngine, l'Integrator,
i Constraints, e i SimulationModule ci lavorano tramite span.

---

## 4. Backend polymorphism

**Decisione: Strategy pattern (classe astratta virtuale).**

```cpp
class ForceEngineBackend {
public:
    virtual ~ForceEngineBackend() = default;
    virtual EnergyReport compute(
        SystemData& sysdata,
        const Cell& cell,
        const std::vector<std::unique_ptr<ForceComponent>>& components
    ) = 0;
};
```

### Backend concreti previsti

- `CPUBackend` — implementazione di riferimento, singolo thread (poi OpenMP)
- `CUDABackend` — NVIDIA GPU, kernel CUDA
- `MetalBackend` — Apple Silicon, kernel Metal
- `OpenCLBackend` — AMD/Intel, kernel OpenCL

### Selezione a init

Il backend viene scelto in base all'hardware disponibile e iniettato nel
ForceEngine prima della prima `compute()`. Il costruttore di ogni backend
può allocare buffer GPU persistenti o compilare kernel.

### Nota su CRTP e std::variant

Vedi Appendice A per le alternative discardate e il loro rationale formativo.

---

## 5. Topologia e parametri risolti

Il `ForceFieldPlugin` produce, per ogni componente, un sottoblocco di
parametri con array SoA già pronti. Non esistono lookup a runtime.

### HarmonicBond

```cpp
struct HarmonicBondParams {
    size_t n;
    std::span<const int>    i, j;
    std::span<const double> k, r0;
};
```
Energia: `E = 0.5 * k * (r - r0)^2`

### HarmonicAngle

```cpp
struct HarmonicAngleParams {
    size_t n;
    std::span<const int>    i, j, k;       // j = atomo centrale
    std::span<const double> k_theta, theta0;
};
```
Energia: `E = 0.5 * k_theta * (theta - theta0)^2`

### PeriodicDihedral

```cpp
struct PeriodicDihedralParams {
    size_t n;
    std::span<const int>    i, j, k, l;
    std::span<const int>    periodicity;
    std::span<const double> k_phi, phi0;
};
```
Energia: `E = k_phi * (1 + cos(n * phi - phi0))`

Nota: i diedri impropers sono inclusi qui con `periodicity = 0` o con un
flag nel tipo.

### LennardJones

```cpp
struct LennardJonesParams {
    size_t n_types;
    std::vector<double> sigma;    // [n_types][n_types] row-major
    std::vector<double> epsilon;  // [n_types][n_types] row-major

    // Combinatorial rules (Lorentz-Berthelot)
    double combine_epsilon(int ti, int tj) const;
    double combine_sigma(int ti, int tj) const;
};
```

Nota: lookup table pre-calcolata. Le regole di combinazione si applicano
una volta dal ForceFieldPlugin, non a runtime.

### CoulombDirect + PME

```cpp
struct CoulombParams {
    // Direct sum
    double cutoff;                        // raggio di troncamento
    double ewald_coefficient;             // alpha per Ewald split

    // PME grid
    struct PMEGrid {
        int nx, ny, nz;
        int spline_order;                 // default 4
        FFTW_Plan fft_forward;
        FFTW_Plan fft_backward;
    } pme;
};
```

### Constraints

```cpp
struct ConstraintsParams {
    size_t n;
    std::span<const int>    i, j;
    std::span<const double> distance_target;
    ConstraintAlgorithm algorithm;        // enum: shake, rattle, settle
};
```

---

## 6. PME + FFTW

### Dipendenza esterna

- CPU: **FFTW 3.x** (o FFTW 3.3+)
- GPU: **cuFFT** (NVIDIA), **clFFT** (AMD via OpenCL)

### Parametri

| Parametro | Default | Descrizione |
|-----------|---------|-------------|
| `cutoff` | 1.2 nm | Raggio troncamento direct sum |
| `ewald_coefficient` | 3.5 / cutoff | Split parameter alpha |
| `grid_spacing` | 0.1 nm | Risoluzione griglia PME |
| `spline_order` | 4 | Interpolazione B-spline |

### Split del potenziale

```
Coulomb_ij = erfc(alpha * r_ij) / r_ij   [direct sum, cutoff]
           + erf(alpha * r_ij) / r_ij     [reciprocal sum, FFT]
           - self term                     [correzione auto-interazione]
```

Ogni termine è calcolato dal rispettivo componente:
- `CoulombDirect`: solo erfc/r nel cutoff (stessa Verlet list di LJ)
- `CoulombPME`: griglia → FFT → B-spline interpolation
- Correzione excluded pairs (1-2, 1-3) intra-molecolari

---

## 7. Precisione

`double` IEEE 754 a 64 bit per tutte le quantità fisiche:
- coordinate e velocità
- forze
- energia e suoi contributi
- parametri di potenziale

Motivazione: la stabilita` a lungo termine di traiettorie MD (milioni di
step) richiede accumulo a 64 bit per evitare derive energetiche.

---

## 8. Error handling

```cpp
enum class ForceError {
    NaN,             // posizioni o forze contengono NaN
    BondExplosion,   // distanza legame > soglia (es. 10x r0)
    Mismatch,        // dimensione span non allineata
    ParmError,       // parametro fisicamente impossibile (sigma <= 0)
};
```

Ogni `compute()` può lanciare `ForceException(ForceError, component_name, context)`.
La eccezione propaga al SimulationEngine che la cattura, scrive un checkpoint
di emergenza e abortisce con messaggio diagnostico completo.

Nel backend GPU si aggiunge:
```cpp
class GPUError : public ForceException {
    std::string kernel_name;
    int device_id;
    cudaError_t cuda_err;  // o equivalente Metal/OpenCL
};
```

---

## 9. Dipendenze software

| Dipendenza | Versione | Uso |
|------------|----------|-----|
| C++ standard | 20 | span, variant, concepts |
| FFTW | 3.x | CPU PME |
| cuFFT | (in CUDA SDK) | GPU PME (NVIDIA) |
| pybind11 | 2.12+ | bridge Python/C++ |

Nessuna dipendenza da Eigen, Boost, o altre librerie di template lineari.

---

## Appendice A — Alternative discardate per Backend polymorphism

### B: CRTP (Curiously Recurring Template Pattern)

```cpp
template <typename Derived>
class ForceEngineBackendBase {
public:
    EnergyReport compute(SystemData& sys, ...) {
        return static_cast<Derived*>(this)->compute_impl(sys, ...);
    }
};

class CPUBackend : public ForceEngineBackendBase<CPUBackend> {
    EnergyReport compute_impl(SystemData& sys, ...) { /* ... */ }
};
```

**Vantaggio:** zero costi di indirezione virtuale — il dispatch è risolto a
compile-time tramite `static_cast<Derived>`.

**Svantaggio:** `template infection`. Tutto ciò che usa il backend come tipo
astratto deve diventare template. Pybind11 non espone classi template senza
esplicita istanziazione. Se si volesse scambiare backend a runtime (es. per
test), serve comunque una classe wrapper virtuale in cima.

### C: Type erasure con std::variant + std::visit

```cpp
using AnyBackend = std::variant<CPUBackend, CUDABackend, MetalBackend>;

struct ComputeVisitor {
    SystemData& sys;
    EnergyReport operator()(CPUBackend& b)   { return b.compute(sys); }
    EnergyReport operator()(CUDABackend& b)  { return b.compute(sys); }
    EnergyReport operator()(MetalBackend& b) { return b.compute(sys); }
};

// dispatch:
EnergyReport e = std::visit(ComputeVisitor{sys}, backend);
```

**Vantaggio:** zero indirezione virtuale (il compilatore genera un jump table
per i casi del variant, simile a uno switch). Idea pulita, moderna, no
polimorfismo dinamico.

**Svantaggio:** il numero di backend è hardcoded nel `variant`. Aggiungere un
nuovo backend significa modificare tutti i visitor (o usare `auto&&` overload
generico). Più verboso dello strategy pattern, e il guadagno in performance
rispetto alla vtable è trascurabile per un metodo che impiega millisecondi
per chiamata.

---

## Appendice B — Riferimenti

- ARCHITECTURE.md — documento architetturale primario del progetto DMD
- GROMACS Reference Manual (2024) — PME, LINCS, Verlet list scheme
- NAMD User's Guide (2023) — orchestrazione componenti forze
- Essmann et al., JCP 103:8577 (1995) — PME formulation
- Ryckaert et al., JCP 23:327 (1977) — RATTLE algorithm
