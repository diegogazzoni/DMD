#include "sim/simulation_engine.h"
#include "sim/simulation_config.h"
#include "core/cell.h"
#include "force/lennard_jones.h"
#include <gtest/gtest.h>
#include <cmath>
#include <filesystem>
#include <vector>
#include <random>

static constexpr double AR_SIGMA = 0.3405;
static constexpr double AR_EPSILON = 0.997;
static constexpr double AR_MASS = 39.948;
static constexpr double KB = 0.008314462618;

static SimulationConfig make_argon_config(double box_nm, size_t n_atoms) {
    SimulationConfig cfg;
    cfg.dt = 0.002;
    cfg.n_steps = 200;
    cfg.box_size = box_nm;
    cfg.use_lj = true;
    cfg.lj = LennardJonesParams{};
    cfg.lj.n_types = 1;
    cfg.lj.sigma = {AR_SIGMA};
    cfg.lj.epsilon = {AR_EPSILON};
    cfg.lj_cutoff = box_nm * 0.49;
    cfg.checkpoint_interval = 0;
    cfg.masses.assign(n_atoms, AR_MASS);
    cfg.charges.assign(n_atoms, 0.0);
    return cfg;
}

static void fcc_lattice(std::vector<double>& x, std::vector<double>& y, std::vector<double>& z,
                        double box_size, int n_cells) {
    double a = box_size / n_cells;
    for (int i = 0; i < n_cells; ++i) {
        for (int j = 0; j < n_cells; ++j) {
            for (int k = 0; k < n_cells; ++k) {
                double bx = i * a, by = j * a, bz = k * a;
                x.push_back(bx);          y.push_back(by);          z.push_back(bz);
                x.push_back(bx + a*0.5);  y.push_back(by + a*0.5);  z.push_back(bz);
                x.push_back(bx + a*0.5);  y.push_back(by);          z.push_back(bz + a*0.5);
                x.push_back(bx);          y.push_back(by + a*0.5);  z.push_back(bz + a*0.5);
            }
        }
    }
}

TEST(CheckpointRestart, TrajectoryReproducibility) {
    std::string ckpt_path = "gate_ckpt.bin";
    int n_cells = 2;
    size_t n_atoms = 4 * n_cells * n_cells * n_cells; // 32
    double density = 1.4;
    double mass_total_g = n_atoms * AR_MASS / 6.02214076e23;
    double box_nm = std::pow(mass_total_g / density, 1.0 / 3.0) * 1e7;

    // initial velocities
    std::mt19937_64 rng(42);
    std::normal_distribution<double> ndist(0.0, 1.0);
    double sigma_v = std::sqrt(KB * 85.0 / AR_MASS);
    std::vector<double> vx_init(n_atoms), vy_init(n_atoms), vz_init(n_atoms);
    for (size_t i = 0; i < n_atoms; ++i) {
        vx_init[i] = ndist(rng) * sigma_v;
        vy_init[i] = ndist(rng) * sigma_v;
        vz_init[i] = ndist(rng) * sigma_v;
    }
    double com_vx = 0, com_vy = 0, com_vz = 0;
    for (size_t i = 0; i < n_atoms; ++i) {
        com_vx += vx_init[i]; com_vy += vy_init[i]; com_vz += vz_init[i];
    }
    com_vx /= n_atoms; com_vy /= n_atoms; com_vz /= n_atoms;
    for (size_t i = 0; i < n_atoms; ++i) {
        vx_init[i] -= com_vx; vy_init[i] -= com_vy; vz_init[i] -= com_vz;
    }

    // positions
    std::vector<double> px, py, pz;
    fcc_lattice(px, py, pz, box_nm, n_cells);

    // ---- Full run of 200 steps ----
    std::vector<double> final_ref;
    {
        auto cfg = make_argon_config(box_nm, n_atoms);
        cfg.pos_x = px; cfg.pos_y = py; cfg.pos_z = pz;
        cfg.vel_x = vx_init; cfg.vel_y = vy_init; cfg.vel_z = vz_init;

        auto engine = build_simulation(cfg);
        engine.run();

        auto& sys = engine.system();
        final_ref.insert(final_ref.end(), sys.pos_x.begin(), sys.pos_x.end());
        final_ref.insert(final_ref.end(), sys.pos_y.begin(), sys.pos_y.end());
        final_ref.insert(final_ref.end(), sys.pos_z.begin(), sys.pos_z.end());
        final_ref.insert(final_ref.end(), sys.vel_x.begin(), sys.vel_x.end());
        final_ref.insert(final_ref.end(), sys.vel_y.begin(), sys.vel_y.end());
        final_ref.insert(final_ref.end(), sys.vel_z.begin(), sys.vel_z.end());
    }

    // ---- Half run: 100 steps, save checkpoint, continue to 200 ----
    std::vector<double> final_restart;
    {
        auto cfg = make_argon_config(box_nm, n_atoms);
        cfg.n_steps = 100;
        cfg.pos_x = px; cfg.pos_y = py; cfg.pos_z = pz;
        cfg.vel_x = vx_init; cfg.vel_y = vy_init; cfg.vel_z = vz_init;

        auto engine = build_simulation(cfg);
        engine.run();
        EXPECT_EQ(engine.current_step(), 100);
        engine.save_checkpoint(ckpt_path);
    }

    // fresh engine, load checkpoint, run remaining 100 steps
    {
        auto cfg = make_argon_config(box_nm, n_atoms);
        cfg.n_steps = 100;
        cfg.pos_x = px; cfg.pos_y = py; cfg.pos_z = pz;
        cfg.vel_x = vx_init; cfg.vel_y = vy_init; cfg.vel_z = vz_init;

        auto engine = build_simulation(cfg);
        engine.load_checkpoint(ckpt_path);
        EXPECT_EQ(engine.current_step(), 100);
        engine.run();
        EXPECT_EQ(engine.current_step(), 200);

        auto& sys = engine.system();
        final_restart.insert(final_restart.end(), sys.pos_x.begin(), sys.pos_x.end());
        final_restart.insert(final_restart.end(), sys.pos_y.begin(), sys.pos_y.end());
        final_restart.insert(final_restart.end(), sys.pos_z.begin(), sys.pos_z.end());
        final_restart.insert(final_restart.end(), sys.vel_x.begin(), sys.vel_x.end());
        final_restart.insert(final_restart.end(), sys.vel_y.begin(), sys.vel_y.end());
        final_restart.insert(final_restart.end(), sys.vel_z.begin(), sys.vel_z.end());
    }

    // ---- Compare ----
    ASSERT_EQ(final_ref.size(), final_restart.size());
    for (size_t i = 0; i < final_ref.size(); ++i) {
        EXPECT_DOUBLE_EQ(final_ref[i], final_restart[i]) << "Mismatch at " << i;
    }

    std::filesystem::remove(ckpt_path);
}
