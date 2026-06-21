# DMD — C4 Level 3: Component Diagram

## Simulation Step Loop

```mermaid
graph TB
    subgraph Step["SimulationEngine::step() — Single MD step"]
        direction TB
        ZF1["1. Zero forces"]
        FC1["2. Compute forces (LJ + Coulomb + Bonded)"]
        HK1["3. Half-kick (v += F/m · dt/2)"]
        ADV["4. Advance positions (x += v · dt)"]
        SHK["5. SHAKE constraint correction"]
        ZF2["6. Zero forces"]
        FC2["7. Compute forces (same components)"]
        HK2["8. Half-kick (v += F/m · dt/2)"]
        RTL["9. RATTLE velocity correction"]
        THM["10. Thermostat apply"]
        BAR["11. Barostat apply"]
        INC["12. step++, time += dt"]
        
        ZF1 --> FC1 --> HK1 --> ADV --> SHK --> ZF2 --> FC2 --> HK2 --> RTL --> THM --> BAR --> INC
    end

    subgraph Forces["Force Components"]
        LJ["LennardJones
            E = 4ε[(σ/r)¹² - (σ/r)⁶]
            Verlet neighbor list"]
        HB["HarmonicBond
            E = ½k(r - r₀)²"]
        HA["HarmonicAngle
            E = ½k_θ(θ - θ₀)²"]
        PD["PeriodicDihedral
            E = k_φ[1+cos(nφ-φ₀)]"]
        CD["CoulombDirect
            E = k_e·q_iq_j·erfc(αr)/r"]
        PME["CoulombPME
            B-spline → FFT → reciprocal"]
        CE["CoulombExclusion
            Correction for 1-2/1-3/1-4 pairs"]
    end

    subgraph Ensemble["Ensemble Controllers"]
        NH["NoseHooverThermostat
            ζ̇ = (T - T₀)/Q · N_f · k_B"]
        BR["BerendsenThermostat
            v ← v · √(1 + dt/τ(T/T₀ - 1))"]
        AN["AndersenThermostat
            Stochastic collisions at rate ν"]
        
        BBar["BerendsenBarostat
            x ← x · (1 - dt/τ·(P₀-P)/β)"]
        ABar["AndersenBarostat
            Extended piston ϵ̇, W"]
    end

    subgraph I_O["I/O"]
        Traj["H5MDWriter
            /particles/atoms/{position,velocity,box}"]
        CKPT["Checkpoint
            Binary: magic + system + Verlet state"]
        Progress["SimulationStatus
            GROMACS-style progress table"]
    end

    FC1 -.- Forces
    FC2 -.- Forces
    THM -.- Ensemble
    BAR -.- Ensemble
    Step -.- I_O
```

## Python Data Flow

```mermaid
graph LR
    subgraph Input["Input Preparation"]
        PDB["from_pdb()
            PDB → positions + box"]
        PSF["from_psf()
            PSF → topology + exclusions"]
        FF["load_forcefield()
            CHARMM/AMBER → FFParams"]
    end
    
    subgraph Merge["Force Field Merge"]
        MERGE["merge_ff()
            Match PSF types → FF params"]
        CONSTRS["generate_constraints()
            PSF angle → SHAKE pairs"]
    end
    
    subgraph System["System Construction"]
        SB["SystemBuilder
            dict → SimulationConfig"]
        APP["apply_config()
            config.json → SimulationConfig"]
    end
    
    subgraph Run["Execution"]
        RUN["run()
            Engine(cfg) + progress"]
        ENG["EngineWrapper
            C++ engine via pybind11"]
    end
    
    subgraph Output["Output"]
        H5MD["trajectory.h5
            H5MD standard"]
        ENERGY["energy.h5
            Step-by-step energy"]
        CKPT["checkpoint.bin
            Binary restart"]
    end

    PDB --> SB
    PSF --> MERGE
    PSF --> CONSTRS
    FF --> MERGE
    FF --> CONSTRS
    MERGE --> SB
    CONSTRS --> SB
    SB --> RUN
    APP --> RUN
    RUN --> ENG
    ENG --> H5MD
    ENG --> ENERGY
    ENG --> CKPT
```

## Data Structures

| Struct | Layout | Note |
|--------|--------|------|
| **SystemData** | SoA: `pos_x/y/z`, `vel_x/y/z`, `forces_x/y/z`, `masses`, `charges`, `atom_types` | Contiguous arrays, `std::span` accessors |
| **SimulationConfig** | Flat struct: scalars + `std::vector<T>` | Populated by Python, consumed by `build_simulation()` |
| **SimulationEngine::Config** | Subset: dt, n_steps, trajectory/path settings | Used by engine at runtime |
