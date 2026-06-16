---
name: dmd-workflow
description: >
  Ciclo di sviluppo DMD: implementa, testa, commit, aggiorna PROGRESS.md.
  Ogni feature segue: design, impl, test, build, ctest, commit, report, next.
  Usa per ogni feature aggiunta al motore MD.
metadata:
  project: dmd
  phase: "1"
---

## Ciclo di sviluppo

### Sequenza obbligatoria per ogni feature

```
 1. DESIGN  → definisci interfaccia/struttura negli header
 2. IMPL    → implementa nel .cpp
 3. TEST    → scrivi test in tests/unit/<modulo>/test_*.cpp
 4. BUILD   → cmake --build && ctest --output-on-failure
 5. BUGFIX  → se test falliscono, correggi e ripeti BUILD
 6. COMMIT  → git add -A && git commit -m "feat(<modulo>): messaggio"
 7. REPORT  → aggiorna PROGRESS.md
 8. NEXT    → passa alla feature successiva
```

### Regole per i test
- Ogni file sorgente ha un corrispondente test unitario
- Copri caso nominale + edge case (zero, negativo, NaN)
- Usa `EXPECT_DOUBLE_EQ` per confronti floating-point
- Mock di ciò che sta sopra il layer testato (`using ::testing::StrictMock`)

### Regole per il commit
- Messaggio in inglese, formato conventional commit
- Esempi:
  - `feat(core): add Cell with 4 types and minimum_image`
  - `fix(force): correct sign in harmonic bond force`
  - `test(force): add three-body angle test`
- Non committare mai file binari, .DS_Store, o test che non passano

### Regole per PROGRESS.md
- Usa formato cronologico (data come titolo)
- Per ogni entry: file creati/modificati, test passati, gate raggiunto
- Aggiungi note tecniche o decisioni cambiate
- Template entry:

```markdown
### 2026-06-16 — Nome feature

- File: `src/modulo/file.h`, `src/modulo/file.cpp`, `tests/unit/modulo/test_file.cpp`
- Test: N/M passati
- Gate: —
- Note: [eventuali osservazioni]
```
