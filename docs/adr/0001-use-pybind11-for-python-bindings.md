# ADR 0001: Use pybind11 for Python bindings

**Status:** Accepted

**Context:** DMD's C++ simulation core needs to be callable from Python. Options included pybind11, Cython, ctypes, and SWIG.

**Decision:** Use pybind11 because:
- Minimal boilerplate (modern C++ reflection, header-only)
- First-class NumPy array support (zero-copy for forces, positions, etc.)
- Active maintenance, CMake FetchContent integration
- Used by major projects (e.g., OpenAI Gym, MFEM)

**Consequences:**
- C++ headers must be pybind11-compatible (no templates in exposed APIs)
- Build requires a C++17 compiler with RTTI
- Each new SimulationConfig field needs a `.def_readwrite()` line
