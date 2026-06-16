# ARCHITECTURE.md

Documento di riferimento primario per il motore di dinamica molecolare modulare. Qualsiasi agente AI che operi su questo repository deve leggere questo file prima di scrivere o modificare codice. In caso di conflitto tra questo documento e un'implementazione esistente, il documento ha precedenza: l'implementazione va corretta, non il documento.

Le interfacce descritte qui sono definite formalmente in `md_engine_core/interfaces.py` e verificate da `md_engine_core/test_interfaces.py`. Lo schema di validazione strutturale del file di sistema è in `msys.schema.json`. Questi tre file, insieme al presente documento, costituiscono l'intera Fase 0 del progetto.

---

## 1. Architettura a layer

Il motore è organizzato in layer sovrapposti. Ogni layer conosce esclusivamente le interfacce del layer immediatamente sottostante.

```
API Layer            Python SDK, REST API, CLI
        |
Simulation Engine     coordina il ciclo MD: Integrator, EnsembleController
        |
Module Registry       ForceFieldPlugin (init) e SimulationModule (hook)
        |
Force Engine          bonded, non-bonded, long-range, constraint; backend CPU/CUDA/Metal/OpenCL
        |
System Layer          parsing JSON, compilazione binaria interna, Cell, I/O (H5MD, checkpoint)
```

**Regola di dipendenza (invariante I-2):** un layer non importa mai un'implementazione concreta di un layer superiore, e non salta layer intermedi per "comodità". Il Force Engine non conosce l'Integrator. Il Module Registry non conosce l'API Layer. Se un agente si trova a importare qualcosa che viola questa direzione, l'architettura è sbagliata, non il caso d'uso.

---

## 2. Contratti tra layer

| Interfaccia | Input | Output | Definita in |
|---|---|---|---|
| `Cell` | — | wrap/minimum_image/volume/scale, matrice 3x3 | `interfaces.py` |
| `ForceEngine` | posizioni, topologia, cella | forze, energia potenziale | `interfaces.py` |
| `Integrator` | pos, vel, forze, masse, dt | pos, vel aggiornate | `interfaces.py` |
| `ThermostatAlgorithm` / `BarostatAlgorithm` | vel/cella, target T/P, dt | vel/cella aggiornate + stato esteso | `interfaces.py` |
| `EnsembleController` | percorso preset, parametri | istanza di algoritmo caricata | `interfaces.py` |
| `ForceFieldPlugin` | `RawSystem` | `ParametrizedSystem` | `interfaces.py` |
| `SimulationModule.hook` | `SimContext` (letto) | modifiche tramite metodi espliciti | `interfaces.py` |
| `SystemLayer.read_system` | path JSON | `RawSystem` (già validato) | `interfaces.py` |
| `SystemLayer` checkpoint/trajectory | `CheckpointState` / `SimContext` | file checkpoint JSON / frame H5MD | `interfaces.py` |

Ogni interfaccia è verificata in `test_interfaces.py` con un'implementazione dummy minimale, per garantire che il contratto sia effettivamente implementabile prima di passare a un'implementazione di produzione.

---

## 3. Ecosistema file

L'utente produce **solo JSON**. Non esiste alcun formato binario visibile all'utente, in input o in output del flusso utente.

```
system.json  →  [MD Engine]  →  trajectory.h5md
                              →  checkpoint.json (periodico + sempre all'uscita)
```

Il file di sistema (`system.json`) è unificato in tre sezioni — `system`, `force_field`, `simulation` — secondo lo schema strutturale in `msys.schema.json`. Il force field **non è mai inline**: è sempre un riferimento a un plugin risolto dal Module Registry (`force_field.plugin` + `force_field.version`). Questo evita la duplicazione di force field standard (Amber, CHARMM, OpenFF) su migliaia di sistemi diversi.

L'output di traiettoria è **H5MD nativo** (HDF5). Nessun formato proprietario: i converter verso DCD, XTC o altri formati GROMACS-compatibili sono strumenti separati, non parte del core, perché MDAnalysis e l'ecosistema esistente leggono H5MD direttamente.

---

## 4. Formato binario interno

Il motore compila `system.json` in una rappresentazione binaria interna al momento dell'avvio. Questo binario **non è mai scritto su disco come artefatto utente**: è un dettaglio implementativo del System Layer, usato per allocare le strutture dati in modo efficiente prima di entrare nel ciclo MD.

```
[ Header: MAGIC · VERSION · offset alle sezioni ]
[ System section: box matrix, cell type, n_atoms, molecole ]
[ FF section: tipi atomici, LJ, bonded, cariche ]
[ Positions section: coordinate, velocità (opzionali / restart) ]
```

L'header contiene gli offset diretti a ogni sezione, per permettere seek senza lettura sequenziale (es. un restart che necessita solo della sezione posizioni). Lo stesso formato è condiviso tra avvio nuovo e restart: cambia solo la provenienza della sezione FF (output del `ForceFieldPlugin` oppure deserializzata dal checkpoint).

---

## 5. Cella di simulazione

Rappresentazione interna sempre come matrice 3x3. Quattro tipi concreti, stessa interfaccia `Cell`:

- **Ortorombica** — caso generale, matrice diagonale.
- **Esagonale** — membrane, canali ionici. `minimum_image()` richiede 7 immagini candidate (la cella primaria non è convessa in proiezione XY); usare l'algoritmo ortorombico su una cella esagonale produce errori fisici silenziosi.
- **Triclina** — caso generale, matrice arbitraria.
- **Dodecaedro rombico** — proteine globulari solvatate, packing efficiente.

Il `ForceEngine` e l'`Integrator` dipendono solo dall'interfaccia `Cell`, mai dal tipo concreto. Il checkpoint serializza sempre la matrice 3x3 completa, indipendentemente dal tipo, così un restart è riproducibile anche se il tipo dichiarato cambia.

---

## 6. Ensemble: termostati e barostati come preset

Termostati e barostati non sono hardcoded nell'`EnsembleController` (invariante I-3): sono algoritmi (`ThermostatAlgorithm`, `BarostatAlgorithm`) caricati da file preset esterni.

```
thermostats/velocity_rescaling.preset
thermostats/berendsen.preset
thermostats/nose_hoover.preset
barostats/berendsen.preset
barostats/parrinello_rahman.preset
```

Ogni algoritmo dichiara le proprie `state_variables` (es. Nosé-Hoover: `xi`, `xi_dot`). Queste sono serializzate automaticamente nel checkpoint tramite `get_state()` / `set_state()` (invariante I-5), senza che l'`EnsembleController` debba conoscerne la struttura interna.

Aggiungere un nuovo algoritmo significa scrivere un nuovo preset e la classe che implementa `ThermostatAlgorithm` o `BarostatAlgorithm`, senza toccare l'architettura esistente.

---

## 7. Module Registry: due tipi di plugin distinti

Il Module Registry gestisce due categorie di plugin con interfacce e momenti di esecuzione diversi. Non vanno confusi.

**`ForceFieldPlugin`** — eseguito una sola volta, in fase di inizializzazione, prima che il ciclo MD inizi. Espone `parametrize(system) -> ParametrizedSystem`. Deve fallire esplicitamente (eccezione) se la topologia non è completamente coperta — atomo non riconosciuto, legame non parametrizzato — mai con un default silenzioso.

**`SimulationModule`** — si aggancia agli hook del ciclo MD (`on_init`, `pre_force`, `post_force`, `post_step`, `on_finalize`). Tutti i metodi hanno default no-op: un modulo implementa solo gli hook che gli servono, senza overhead per quelli inutilizzati. Esempi: `ExternalForce` (campo elettrico uniforme), enhanced sampling (umbrella sampling, metadinamica), REMD. L'agente di simulazione reattivo (protocolli di esperimento, vedi roadmap originale) è anch'esso un `SimulationModule`, attualmente in stand-by di progettazione.

```
for step in range(n_steps):
    modules.pre_force()
    forces = force_engine.compute(...)
    modules.post_force()
    positions, velocities = integrator.step(...)
    ensemble_controller.apply(...)
    modules.post_step()
    [write trajectory]
    [write checkpoint]
```

---

## 8. Force Engine e backend GPU

Il Force Engine delega il calcolo a un `ForceEngineBackend` selezionato ad inizializzazione in base all'hardware disponibile:

- **CPU** — implementazione di riferimento.
- **CUDA** — NVIDIA.
- **Metal (MPS)** — Apple Silicon, via `metal-cpp`.
- **OpenCL** — AMD, Intel, portabile.

Linguaggio core: **C++20**. È l'unico linguaggio con binding nativi maturi su tutti e tre gli stack GPU contemporaneamente (CUDA C++ è un'estensione diretta di C++; Metal ha binding C++ ufficiali; OpenCL ha binding C++ standard). Il resto del motore (API, CLI, SDK, orchestrazione) resta in Python, con `pybind11` come bridge. Il loop critico di calcolo forze gira interamente in C++: il garbage collector di Python non lo tocca, perché non vengono creati oggetti Python nel percorso caldo.

---

## 9. Validazione a due passaggi

Il file di sistema attraversa due validazioni sequenziali e bloccanti, mai fuse in un solo passaggio:

```
load JSON
  → JSON Schema validation (msys.schema.json)
      fail → abort, errore di forma
      pass → Pydantic parsing + validatori custom
                fail → abort, errore semantico con riferimento al campo
                pass → RawSystem in memoria
```

**JSON Schema** valida la forma: tipi, campi obbligatori, enum, range numerici, vincoli condizionali tra blocchi (es. NPT richiede `barostat` e `thermostat`; NVE li proibisce esplicitamente). Questo livello non richiede contesto esterno al documento.

**Pydantic con validatori custom** valida la semantica che richiede contesto esterno o cross-referenze: che il plugin force field dichiarato esista nel registry, che gli indici atomici nei legami siano validi, che i parametri di un preset siano coerenti con l'algoritmo dichiarato.

Fondere i due passaggi è un errore: senza la validazione strutturale come gate preliminare, i validatori Pydantic riceverebbero documenti malformati e produrrebbero errori tecnici (`KeyError`) invece di messaggi diagnostici chiari.

---

## 10. Checkpoint autocontenuto

Il checkpoint è scritto in JSON, periodicamente (frequenza configurabile in `simulation.checkpoint.frequency_steps`) e sempre all'uscita, sia graceful sia su errore o segnale.

Contiene la **topologia parametrizzata completa**, non solo coordinate e velocità (invariante I-4). La sezione `provenance` registra quale plugin force field ha generato la topologia, ma è puramente informativa: il motore non la rilegge. Questo significa che un restart funziona anche se il plugin force field originale è stato modificato o rimosso — il checkpoint non dipende da nulla al di fuori di se stesso.

```
Avvio nuovo:    JSON → ForceFieldPlugin.parametrize() → topologia in memoria → MD loop
Restart:        checkpoint.json → topologia in memoria → MD loop (FF plugin non invocato)
```

---

## 11. Invarianti

Queste regole non vanno violate da nessuna implementazione concreta, indipendentemente da quanto un'alternativa sembri più comoda in un caso specifico.

- **I-1** — `SystemLayer` è l'unico componente che tocca l'I/O. `ForceEngine` non scrive mai su disco.
- **I-2** — ogni layer conosce solo il layer immediatamente sottostante (regola di dipendenza, sezione 1).
- **I-3** — termostati e barostati non sono mai hardcoded nell'`EnsembleController`: sempre caricati da preset esterni.
- **I-4** — il checkpoint è sempre autocontenuto: include la topologia parametrizzata, non richiede il `ForceFieldPlugin` originale per il restart.
- **I-5** — le variabili estese di un algoritmo (`state_variables`) sono dichiarate dall'algoritmo stesso e serializzate automaticamente nel checkpoint, senza che il controller ne conosca la struttura.
- **I-6** — i `SimulationModule` modificano lo stato di simulazione solo tramite i metodi espliciti di `SimContext` (es. `add_force_contribution`), mai mutando direttamente gli array esterni.

---

## 12. Cosa non fare

- Non far scrivere la traiettoria o il checkpoint direttamente dal `ForceEngine` "per semplicità": viola I-1.
- Non fondere validazione strutturale e semantica in un solo passaggio (sezione 9).
- Non hardcodare un algoritmo di termostato o barostato nell'`EnsembleController`: viola I-3.
- Non serializzare nel checkpoint solo coordinate e velocità: un checkpoint che dipende dal `ForceFieldPlugin` originale viola I-4.
- Non far mutare a un `SimulationModule` gli array di `SimContext` direttamente: viola I-6.
- Non implementare un nuovo tipo di cella riusando `minimum_image()` di un altro tipo "perché sembra simile" (es. ortorombico su esagonale): produce errori fisici silenziosi, non un'eccezione.
- Non aggiungere parametri force field inline in `system.json`: il force field è sempre un riferimento a un plugin (sezione 3).

---

## 13. Roadmap di sviluppo

```
Fase 0 — Specifiche e interfacce           non delegabile agli agenti (questo documento + schema + ABC + fixture)
Fase 1 — System layer                       parser, binary compiler, Cell, H5MD writer, checkpoint
Fase 2 — Force engine, backend CPU          gate: argon NVE, energia conservata su 10^6 step
Fase 3 — Simulation engine                  gate: run completa e restart producono traiettorie identiche
Fase 4 — Ensemble controller e preset       gate: NVT Maxwell-Boltzmann, NPT densità corretta
Fase 5a — Module registry e FF plugin       parallela a 5b
Fase 5b — Backend GPU (CUDA, Metal, OpenCL) parallela a 5a
Fase 6 — API layer                          Python SDK, CLI, REST API
Fase 7 — Simulation module avanzati         ExternalForce, enhanced sampling, REMD, agente (futuro)
```

Ogni gate è una condizione necessaria, non opzionale, per procedere alla fase successiva.

---

## 14. Riferimenti

- `msys.schema.json` — schema di validazione strutturale del file di sistema.
- `md_engine_core/interfaces.py` — contratti formali tra layer.
- `md_engine_core/test_interfaces.py` — verifica dei contratti con implementazioni dummy.
