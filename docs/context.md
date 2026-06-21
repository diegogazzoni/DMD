# DMD — C4 Level 1: Context Diagram

```mermaid
graph TB
    User["👤 User
        (Simulation Scientist)"]
    
    DMD["DMD
        Molecular Dynamics Engine"]
    
    SystemJSON["📄 system.json
        Topology + FF params"]
    
    ConfigJSON["📄 config.json
        Run parameters"]
    
    CheckpointJSON["📄 checkpoint.json
        Simulation state"]
    
    TrajectoryH5["📀 trajectory.h5
        H5MD trajectory"]
    
    EnergyH5["📀 energy.h5
        Energy timeseries"]
    
    CheckpointBin["💾 checkpoint.bin
        Binary state"]
    
    VMD["🔬 VMD / PyMOL
        Visualization"]
    
    MDAnalysis["📊 MDAnalysis / h5py
        Analysis"]
    
    DepsFFTW["📦 FFTW3
        PME computation"]
    
    DepsHDF5["📦 HDF5
        Trajectory storage"]
    
    DepsPybind11["📦 pybind11
        Python bindings"]
    
    DepsGTest["📦 Google Test
        Unit testing"]

    User -->|"writes"| SystemJSON
    User -->|"writes"| ConfigJSON
    User -->|"runs"| DMD
    
    DMD -->|"reads"| SystemJSON
    DMD -->|"reads"| ConfigJSON
    DMD -->|"reads/writes"| CheckpointJSON
    DMD -->|"writes"| TrajectoryH5
    DMD -->|"writes"| EnergyH5
    DMD -->|"writes"| CheckpointBin
    
    User -->|"analyzes"| TrajectoryH5
    TrajectoryH5 -->|"visualize"| VMD
    TrajectoryH5 -->|"analyze"| MDAnalysis
    
    DMD -.->|"depends on"| DepsFFTW
    DMD -.->|"depends on"| DepsHDF5
    DMD -.->|"depends on"| DepsPybind11
    DMD -.->|"depends on"| DepsGTest

    classDef system fill:#1b77b6,color:#fff,stroke:#0a4c7a
    classDef person fill:#08427b,color:#fff,stroke:#052a50
    classDef input fill:#f0e68c,stroke:#bdb76b
    classDef output fill:#98fb98,stroke:#3cb371
    classDef external fill:#d3d3d3,stroke:#808080
    classDef dep fill:#ffb6c1,stroke:#dc143c

    class DMD system
    class User person
    class SystemJSON,ConfigJSON,CheckpointJSON input
    class TrajectoryH5,EnergyH5,CheckpointBin output
    class VMD,MDAnalysis external
    class DepsFFTW,DepsHDF5,DepsPybind11,DepsGTest dep
```

## Elementi

| Elemento | Tipo | Descrizione |
|----------|------|-------------|
| **Utente** | Persona | Scienziato che prepara input e analizza output |
| **DMD** | Sistema | Motore MD C++/Python con binding pybind11 |
| **system.json** | Input | Topologia, posizioni, parametri FF |
| **config.json** | Input | Parametri runtime (dt, ensemble, output) |
| **checkpoint.json** | Input/Output | Checkpoint testuale per restart |
| **trajectory.h5** | Output | Traiettoria in formato H5MD standard |
| **energy.h5** | Output | Serie temporale energia |
| **checkpoint.bin** | Output | Checkpoint binario bit-exact |
| **VMD/PyMOL** | Esterno | Visualizzazione traiettorie |
| **MDAnalysis/h5py** | Esterno | Analisi dati |
| **FFTW3** | Dipendenza | Trasformata di Fourier per PME |
| **HDF5** | Dipendenza | Libreria HDF5 per storage |
| **pybind11** | Dipendenza | Bridge C++/Python |
| **Google Test** | Dipendenza | Framework test C++ |
