#include "simulation_engine.h"
#include "simulation_config.h"
#include "force/force_engine.h"
#include "force/force_component.h"
#include "force/lennard_jones.h"
#include "checkpoint.h"
#include <fstream>
#include <vector>

SimulationEngine::SimulationEngine(std::unique_ptr<ForceEngine> fe, Cell cell, size_t n_atoms, Config config)
    : fe_(std::move(fe))
    , cell_(cell)
    , sys_(n_atoms)
    , config_(std::move(config))
{
}

SimulationEngine::~SimulationEngine() = default;

void SimulationEngine::step() {
    sys_.step = step_;
    sys_.time = time_;

    sys_.forces_x.assign(sys_.n_atoms, 0.0);
    sys_.forces_y.assign(sys_.n_atoms, 0.0);
    sys_.forces_z.assign(sys_.n_atoms, 0.0);
    sys_.potential_energy = 0.0;

    fe_->compute(sys_, cell_);

    integrator_.half_kick(sys_, config_.dt);
    integrator_.advance(sys_, config_.dt);

    sys_.forces_x.assign(sys_.n_atoms, 0.0);
    sys_.forces_y.assign(sys_.n_atoms, 0.0);
    sys_.forces_z.assign(sys_.n_atoms, 0.0);
    sys_.potential_energy = 0.0;

    fe_->compute(sys_, cell_);

    integrator_.half_kick(sys_, config_.dt);

    ++step_;
    sys_.step = step_;
    time_ += config_.dt;
}

void SimulationEngine::run() {
    for (int i = 0; i < config_.n_steps; ++i) {
        step();

        if (config_.checkpoint_interval > 0 && step_ % config_.checkpoint_interval == 0) {
            save_checkpoint(config_.checkpoint_path);
        }
    }
}

int SimulationEngine::find_lj_component() const {
    for (size_t i = 0; i < fe_->components().size(); ++i) {
        if (fe_->components()[i]->name() == "LennardJones") {
            return static_cast<int>(i);
        }
    }
    return -1;
}

void SimulationEngine::save_checkpoint(const std::string& path) const {
    int lj_state = 0;
    int idx = find_lj_component();
    if (idx >= 0) {
        auto* lj = dynamic_cast<LennardJones*>(fe_->components()[idx].get());
        if (lj) lj_state = lj->step_since_rebuild();
    }
    checkpoint::save(path, sys_, lj_state);
}

void SimulationEngine::load_checkpoint(const std::string& path) {
    checkpoint::load(path, sys_, lj_step_since_rebuild_);
    step_ = static_cast<int>(sys_.step);
    time_ = sys_.time;

    int idx = find_lj_component();
    if (idx >= 0) {
        auto* lj = dynamic_cast<LennardJones*>(fe_->components()[idx].get());
        if (lj) lj->set_step_since_rebuild(lj_step_since_rebuild_);
    }
}

// ---- build_simulation ----

std::unique_ptr<ForceEngine> build_force_engine(const SimulationConfig& cfg) {
    auto fe = std::make_unique<ForceEngine>();

    if (cfg.use_lj) {
        auto lj = std::make_unique<LennardJones>(cfg.lj, cfg.lj_cutoff);
        fe->add_component(std::move(lj));
    }

    for (auto& b : cfg.bonds) {
        fe->add_component(std::make_unique<HarmonicBond>(b));
    }
    for (auto& a : cfg.angles) {
        fe->add_component(std::make_unique<HarmonicAngle>(a));
    }
    for (auto& d : cfg.dihedrals) {
        fe->add_component(std::make_unique<PeriodicDihedral>(d));
    }

    return fe;
}

SimulationEngine build_simulation(const SimulationConfig& cfg) {
    auto fe = build_force_engine(cfg);

    auto* lj_ptr = [&]() -> LennardJones* {
        for (auto& c : fe->components()) {
            if (auto* lj = dynamic_cast<LennardJones*>(c.get())) return lj;
        }
        return nullptr;
    }();

    Cell cell(cfg.box_size, cfg.box_size, cfg.box_size);

    size_t n = cfg.pos_x.size();
    SimulationEngine::Config engine_cfg;
    engine_cfg.dt = cfg.dt;
    engine_cfg.n_steps = cfg.n_steps;
    engine_cfg.trajectory_interval = cfg.trajectory_interval;
    engine_cfg.checkpoint_interval = cfg.checkpoint_interval;
    engine_cfg.checkpoint_path = cfg.checkpoint_path;

    SimulationEngine engine(std::move(fe), cell, n, engine_cfg);
    auto& sys = engine.system();

    sys.masses.assign(cfg.masses.begin(), cfg.masses.end());
    sys.charges.assign(cfg.charges.begin(), cfg.charges.end());
    sys.pos_x.assign(cfg.pos_x.begin(), cfg.pos_x.end());
    sys.pos_y.assign(cfg.pos_y.begin(), cfg.pos_y.end());
    sys.pos_z.assign(cfg.pos_z.begin(), cfg.pos_z.end());
    sys.vel_x.assign(cfg.vel_x.begin(), cfg.vel_x.end());
    sys.vel_y.assign(cfg.vel_y.begin(), cfg.vel_y.end());
    sys.vel_z.assign(cfg.vel_z.begin(), cfg.vel_z.end());

    std::vector<int> types(n, 0);
    sys.atom_types = types;

    if (lj_ptr) {
        lj_ptr->set_atom_types(sys.atom_types);
    }

    return engine;
}
