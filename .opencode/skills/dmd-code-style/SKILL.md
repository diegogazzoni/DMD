---
name: dmd-code-style
description: >
  Regole di stile C++ per il motore DMD. Priorità: semplicità, eleganza,
  performance. Crea funzioni solo per codice ripetuto 2+ volte. Stessa
  regola per oggetti. Priorità: 1) velocità 2) memoria. Usa quando scrivi
  o modifichi codice C++ del motore MD.
metadata:
  project: dmd
  phase: "1"
---

## Stile di codice

### Regole cardinali
1. **Niente funzioni per operazioni eseguite una sola volta.** Se un blocco
   è usato in un solo punto, resta inline. Refattorizza in funzione solo
   quando lo stesso blocco compare in 2+ posti.
2. **Stessa regola per classi/struct.** Nessuna gerarchia inutile. Preferisci
   struct semplici a classi con metodi se l'oggetto è solo contenitore di dati.
3. **Performance prima di tutto.** L'allocazione dinamica è un warning implicito.
   Usa stack allocation quando possibile.
4. **Memoria contigua.** Preferisci `std::vector` / `std::span` / array C-style
   su linked list o strutture node-based.

### Convenzioni
- Nomi: `snake_case` per funzioni e variabili, `PascalCase` per classi/enum
- Header: `.h` (no `.hpp`)
- `const` su tutti i parametri che non mutano
- `noexcept` su funzioni che non lanciano eccezioni
- Preferisci `double` a `float` per quantità fisiche
- Niente `auto` per tipi numerici (usa `double`, `int`, esplicito)

### Anti-pattern
- ❌ Funzioni "wrapper" di 2 righe
- ❌ Classi con un solo metodo più lungo della classe stessa
- ❌ `shared_ptr` per dati con ownership singola
- ❌ `dynamic_cast` nel loop caldo
- ❌ Template metaprogramming non necessario

### Cosa fare invece
- Usa `std::span` per passare array (porta size inclusa)
- Usa `unique_ptr` per ownership esclusiva
- Usa struct aggregato per parametri raggruppati
- Inline lambdas per logica usa-e-getta
