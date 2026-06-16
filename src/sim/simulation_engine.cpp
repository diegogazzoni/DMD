#include "simulation_engine.h"
#include "simulation_config.h"
#include "core/constants.h"
#include "force/force_engine.h"
#include "force/force_component.h"
#include "force/lennard_jones.h"
#include "force/coulomb_direct.h"
#include "force/coulomb_pme.h"
#include "thermostat/thermostat.h"
#include "thermostat/berendsen.h"
#include "thermostat/nose_hoover.h"
#include "thermostat/andersen.h"
#include "barostat/barostat.h"
#include "barostat/berendsen_barostat.h"
#include "barostat/andersen_barostat.h"
#include "checkpoint.h"
#include "trajectory/h5md_writer.h"
#include <fstream>
#include <random>
#include <vector>

SimulationEngine::SimulationEngine(
    std::unique_ptr<ForceEngine> fe,
    Cell cell,
    size_t n_atoms,
    Config config,
    std::unique_ptr<Thermostat> thermostat,
    std::unique_ptr<Barostat> barostat
)
    : fe_(std::move(fe))
    , thermostat_(std::move(thermostat))
    , barostat_(std::move(barostat))
    , cell_(cell)
    , sys_(n_atoms)
    , config_(std::move(config))
{
}

SimulationEngine::~SimulationEngine() = default;

void SimulationEngine::set_thermostat(std::unique_ptr<Thermostat> t) {
    thermostat_ = std::move(t);
}

void SimulationEngine::set_barostat(std::unique_ptr<Barostat> b) {
    barostat_ = std::move(b);
}

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

    if (thermostat_) thermostat_->apply(sys_, config_.dt);
    if (barostat_) barostat_->apply(sys_, cell_, config_.dt);

    ++step_;
    sys_.step = step_;
    time_ += config_.dt;
}

void SimulationEngine::run() {
    if (config_.trajectory_interval > 0 && !trajectory_writer_) {
        trajectory_writer_ = std::make_unique<H5MDWriter>(
            config_.trajectory_path, sys_.n_atoms,
            /*write_vel=*/true);
    }

    for (int i = 0; i < config_.n_steps; ++i) {
        step();

        if (config_.trajectory_interval > 0 && step_ % config_.trajectory_interval == 0) {
            save_frame();
        }

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

void SimulationEngine::save_frame() {
    if (trajectory_writer_) {
        trajectory_writer_->write_frame(sys_, cell_);
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

    if (cfg.coulomb_type == "direct" || cfg.coulomb_type == "cutoff") {
        CoulombDirectParams p;
        p.cutoff = cfg.coulomb_cutoff;
        p.ewald_coefficient = cfg.ewald_coeff;
        fe->add_component(std::make_unique<CoulombDirect>(p));
    } else if (cfg.coulomb_type == "pme") {
        PMEParams p;
        p.cutoff = cfg.coulomb_cutoff;
        p.ewald_coefficient = cfg.ewald_coeff;
        p.spline_order = cfg.pme_order;
        p.nx = std::max(8, static_cast<int>(cfg.box_size / cfg.pme_grid_spacing + 0.5));
        p.ny = p.nx;
        p.nz = p.nx;
        fe->add_component(std::make_unique<CoulombPME>(p));
    }

    return fe;
}

static std::unique_ptr<Thermostat> make_thermostat(const SimulationConfig& cfg, size_t n_atoms) {
    auto& t = cfg.thermostat;
    if (t.type == "none") return nullptr;
    if (t.type == "berendsen") {
        BerendsenParams p;
        p.temperature = t.temperature;
        p.tau = t.tau;
        return std::make_unique<BerendsenThermostat>(p);
    }
    if (t.type == "nose-hoover") {
        NoseHooverParams p;
        p.temperature = t.temperature;
        p.tau = t.tau;
        return std::make_unique<NoseHooverThermostat>(p, n_atoms);
    }
    if (t.type == "andersen") {
        AndersenParams p;
        p.temperature = t.temperature;
        p.frequency = t.frequency;
        return std::make_unique<AndersenThermostat>(p);
    }
    throw std::runtime_error("Unknown thermostat type: " + t.type);
}

static std::unique_ptr<Barostat> make_barostat(const SimulationConfig& cfg, size_t n_atoms) {
    auto& b = cfg.barostat;
    if (b.type == "none") return nullptr;
    if (b.type == "berendsen") {
        BerendsenBarostatParams p;
        p.pressure = b.pressure;
        p.tau = b.tau;
        p.compressibility = b.compressibility;
        return std::make_unique<BerendsenBarostat>(p);
    }
    if (b.type == "andersen") {
        AndersenBarostatParams p;
        p.pressure = b.pressure;
        p.tau = b.tau;
        return std::make_unique<AndersenBarostat>(p, n_atoms);
    }
    throw std::runtime_error("Unknown barostat type: " + b.type);
}

SimulationEngine build_simulation(const SimulationConfig& cfg) {
    auto fe = build_force_engine(cfg);

    auto* lj_ptr = [&]() -> LennardJones* {
        for (auto& c : fe->components()) {
            if (auto* lj = dynamic_cast<LennardJones*>(c.get())) return lj;
        }
        return nullptr;
    }();

    // Set charges on Coulomb components before moving fe
    auto* coulomb_direct = [&]() -> CoulombDirect* {
        for (auto& c : fe->components()) {
            if (auto* cd = dynamic_cast<CoulombDirect*>(c.get())) return cd;
        }
        return nullptr;
    }();
    auto* coulomb_pme = [&]() -> CoulombPME* {
        for (auto& c : fe->components()) {
            if (auto* cp = dynamic_cast<CoulombPME*>(c.get())) return cp;
        }
        return nullptr;
    }();

    Cell cell(cfg.box_size, cfg.box_size, cfg.box_size);

    size_t n = cfg.pos_x.size();
    SimulationEngine::Config engine_cfg;
    engine_cfg.dt = cfg.dt;
    engine_cfg.n_steps = cfg.n_steps;
    engine_cfg.trajectory_interval = cfg.nstxout;
    engine_cfg.checkpoint_interval = cfg.checkpoint_interval;
    engine_cfg.trajectory_path = cfg.trajectory_path;
    engine_cfg.checkpoint_path = cfg.checkpoint_path;

    auto thermo = make_thermostat(cfg, n);
    auto baro = make_barostat(cfg, n);

    SimulationEngine engine(std::move(fe), cell, n, engine_cfg,
                            std::move(thermo), std::move(baro));
    auto& sys = engine.system();

    sys.masses.assign(cfg.masses.begin(), cfg.masses.end());
    sys.charges.assign(cfg.charges.begin(), cfg.charges.end());
    sys.pos_x.assign(cfg.pos_x.begin(), cfg.pos_x.end());
    sys.pos_y.assign(cfg.pos_y.begin(), cfg.pos_y.end());
    sys.pos_z.assign(cfg.pos_z.begin(), cfg.pos_z.end());
    sys.vel_x.assign(cfg.vel_x.begin(), cfg.vel_x.end());
    sys.vel_y.assign(cfg.vel_y.begin(), cfg.vel_y.end());
    sys.vel_z.assign(cfg.vel_z.begin(), cfg.vel_z.end());

    sys.atom_types = cfg.atom_types;
    if (sys.atom_types.empty()) {
        sys.atom_types.assign(n, 0);
    }

    if (lj_ptr) {
        lj_ptr->set_atom_types(sys.atom_types);
    }
    if (coulomb_direct) {
        coulomb_direct->set_charges(sys.charges);
    }
    if (coulomb_pme) {
        coulomb_pme->set_charges(sys.charges);
    }

    // Generate velocities if requested
    if (cfg.gen_vel && n > 0) {
        std::mt19937 rng(cfg.seed);
        std::normal_distribution<double> normal(0.0, 1.0);
        double kBT = dmd_constants::KB * cfg.init_temperature;
        for (size_t i = 0; i < n; ++i) {
            double sigma = std::sqrt(kBT / sys.masses[i]);
            sys.vel_x[i] = normal(rng) * sigma;
            sys.vel_y[i] = normal(rng) * sigma;
            sys.vel_z[i] = normal(rng) * sigma;
        }
        // Remove center-of-mass momentum
        double com_vx = 0, com_vy = 0, com_vz = 0;
        double total_mass = 0;
        for (size_t i = 0; i < n; ++i) {
            com_vx += sys.masses[i] * sys.vel_x[i];
            com_vy += sys.masses[i] * sys.vel_y[i];
            com_vz += sys.masses[i] * sys.vel_z[i];
            total_mass += sys.masses[i];
        }
        com_vx /= total_mass;
        com_vy /= total_mass;
        com_vz /= total_mass;
        for (size_t i = 0; i < n; ++i) {
            sys.vel_x[i] -= com_vx;
            sys.vel_y[i] -= com_vy;
            sys.vel_z[i] -= com_vz;
        }
        // Rescale to exact T
        if (cfg.init_temperature > 0) {
            double ekin = 0;
            for (size_t i = 0; i < n; ++i) {
                ekin += 0.5 * sys.masses[i] * (
                    sys.vel_x[i] * sys.vel_x[i] +
                    sys.vel_y[i] * sys.vel_y[i] +
                    sys.vel_z[i] * sys.vel_z[i]);
            }
            double T_actual = temperature(ekin, sys.n_atoms);
            double lambda = std::sqrt(cfg.init_temperature / T_actual);
            for (size_t i = 0; i < n; ++i) {
                sys.vel_x[i] *= lambda;
                sys.vel_y[i] *= lambda;
                sys.vel_z[i] *= lambda;
            }
        }
    }

    return engine;
}
