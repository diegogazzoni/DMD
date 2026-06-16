#pragma once

#include "integrate/integrator.h"
#include "core/cell.h"
#include "core/system_data.h"
#include <memory>
#include <string>

class ForceComponent;
class ForceEngine;
class Thermostat;
class Barostat;
class H5MDWriter;

class SimulationEngine {
public:
    struct Config {
        double dt{0.002};
        int n_steps{1000};
        int trajectory_interval{0};
        int checkpoint_interval{500};
        std::string trajectory_path{"trajectory.h5"};
        std::string checkpoint_path{"checkpoint.bin"};
    };

    SimulationEngine(
        std::unique_ptr<ForceEngine> fe,
        Cell cell,
        size_t n_atoms,
        Config config,
        std::unique_ptr<Thermostat> thermostat = nullptr,
        std::unique_ptr<Barostat> barostat = nullptr
    );
    ~SimulationEngine();

    SimulationEngine(const SimulationEngine&) = delete;
    SimulationEngine& operator=(const SimulationEngine&) = delete;
    SimulationEngine(SimulationEngine&&) noexcept = default;
    SimulationEngine& operator=(SimulationEngine&&) noexcept = default;

    void run();
    void step();

    int current_step() const { return step_; }
    double current_time() const { return time_; }
    SystemData& system() { return sys_; }
    const SystemData& system() const { return sys_; }
    ForceEngine& force_engine() { return *fe_; }

    void set_thermostat(std::unique_ptr<Thermostat> t);
    void set_barostat(std::unique_ptr<Barostat> b);

    void save_checkpoint(const std::string& path) const;
    void load_checkpoint(const std::string& path);

private:
    std::unique_ptr<ForceEngine> fe_;
    std::unique_ptr<Thermostat> thermostat_;
    std::unique_ptr<Barostat> barostat_;
    std::unique_ptr<H5MDWriter> trajectory_writer_;
    Integrator integrator_;
    Cell cell_;
    SystemData sys_;
    Config config_;
    int step_{0};
    double time_{0.0};
    int lj_step_since_rebuild_{0};

    void save_frame();
    int find_lj_component() const;
};
