# ADR 0002: Use PME for long-range electrostatics

**Status:** Accepted

**Context:** Coulomb interactions in periodic systems require handling of long-range contributions. Options included Ewald summation, Particle Mesh Ewald (PME), Particle-Particle Particle-Mesh (P3M), and simple cutoff (reaction field).

**Decision:** Use PME (via FFTW) because:
- O(N log N) scaling vs O(N^2) for direct Ewald
- Industry standard for biomolecular simulation (GROMACS, AMBER, CHARMM)
- Compatible with exclusion corrections for 1-2/1-3/1-4 pairs
- FFTW is widely available and optimized for Apple Silicon

**Consequences:**
- FFTW3 is a hard dependency (Homebrew: `brew install fftw`)
- Grid size and spline order are user-configurable parameters
- Exclusion correction component (CoulombExclusion) required for bonded pairs
