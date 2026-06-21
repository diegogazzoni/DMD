# ADR 0004: Use SHAKE constraints for rigid water models

**Status:** Accepted (implemented 2026-06-18)

**Context:** TIP3P water requires fixed O-H bond lengths and a fixed H-O-H angle. Options included harmonic bonds+angles (flexible) and SHAKE/SETTLE constraints (rigid).

**Decision:** Use SHAKE (position) + RATTLE (velocity) constraints, explicitly enabled via `config.json`:
- Three constraints per water molecule: 2x O-H (0.09572 nm) + 1x H-H (0.15139 nm)
- Harmonic bonds and angles remain available for protein/ligand flexible terms
- Config explicit: `"constraints": {"type": "shake"}` — no hidden auto-detection

**Why not SETTLE:** SETTLE is O(1) and analytical but works only for 3-site water. SHAKE is iterative but general — supports any rigid fragment (water, ammonia, CH₄, protein backbone).

**Consequences:**
- Timestep can increase from 0.5 fs (angle-limited) to 2 fs (SHAKE-limited)
- No bond/angle energy contribution (pure LJ + Coulomb, as TIP3P is calibrated)
- `generate_constraints()` in `merger.py` produces 3 constraints per water from PSF+FF
- SHAKE iteration overhead is negligible vs PME (~2-3 iterations per step)
