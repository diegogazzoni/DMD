#include "sim/simulation_engine.h"
#include "sim/simulation_config.h"
#include "force/lennard_jones.h"
#include <gtest/gtest.h>
#include <cmath>
#include <filesystem>

TEST(SimulationEngineTest, BasicNVELJ) {
    SimulationConfig cfg;
    cfg.dt = 0.002;
    cfg.n_steps = 10;
    cfg.box_size = 3.0;
    cfg.lj_cutoff = 1.2;
    cfg.use_lj = true;
    cfg.lj = LennardJonesParams{};
    cfg.lj.n_types = 1;
    cfg.lj.sigma = {0.3405};
    cfg.lj.epsilon = {0.997};

    int n = 8;
    cfg.masses.assign(n, 39.948);
    cfg.charges.assign(n, 0.0);
    for (int i = 0; i < n; ++i) {
        cfg.pos_x.push_back(0.5 + static_cast<double>(i) * 0.3);
        cfg.pos_y.push_back(0.5 + static_cast<double>(i) * 0.3);
        cfg.pos_z.push_back(0.5 + static_cast<double>(i) * 0.3);
        cfg.vel_x.push_back(0.0);
        cfg.vel_y.push_back(0.0);
        cfg.vel_z.push_back(0.0);
    }

    auto engine = build_simulation(cfg);
    EXPECT_EQ(engine.current_step(), 0);
    EXPECT_EQ(engine.system().n_atoms, 8);

    engine.run();
    EXPECT_EQ(engine.current_step(), 10);
}

TEST(SimulationEngineTest, KineticEnergyPositive) {
    SimulationConfig cfg;
    cfg.dt = 0.001;
    cfg.n_steps = 5;
    cfg.box_size = 4.0;
    cfg.use_lj = true;
    cfg.lj = LennardJonesParams{};
    cfg.lj.n_types = 1;
    cfg.lj.sigma = {0.3405};
    cfg.lj.epsilon = {0.997};
    cfg.lj_cutoff = 2.0;

    int n = 4;
    cfg.masses.assign(n, 1.0);
    cfg.charges.assign(n, 0.0);
    for (int i = 0; i < n; ++i) {
        cfg.pos_x.push_back(1.0 + i);
        cfg.pos_y.push_back(1.0);
        cfg.pos_z.push_back(1.0);
        cfg.vel_x.push_back(0.5);
        cfg.vel_y.push_back(0.0);
        cfg.vel_z.push_back(0.0);
    }

    auto engine = build_simulation(cfg);
    engine.run();

    auto& sys = engine.system();
    double ekin = 0;
    for (size_t i = 0; i < sys.n_atoms; ++i) {
        ekin += 0.5 * sys.masses[i] * (
            sys.vel_x[i]*sys.vel_x[i] +
            sys.vel_y[i]*sys.vel_y[i] +
            sys.vel_z[i]*sys.vel_z[i]);
    }
    EXPECT_GT(ekin, 0);
}

TEST(SimulationEngineTest, SaveAndLoadCheckpoint) {
    std::string ckpt_path = "test_ckpt.bin";

    {
        SimulationConfig cfg;
        cfg.dt = 0.002;
        cfg.n_steps = 5;
        cfg.box_size = 3.0;
        cfg.use_lj = true;
        cfg.lj = LennardJonesParams{};
        cfg.lj.n_types = 1;
        cfg.lj.sigma = {0.3405};
        cfg.lj.epsilon = {0.997};
        cfg.lj_cutoff = 1.2;
        cfg.checkpoint_interval = 5;
        cfg.checkpoint_path = ckpt_path;

        int n = 4;
        cfg.masses.assign(n, 39.948);
        cfg.charges.assign(n, 0.0);
        for (int i = 0; i < n; ++i) {
            cfg.pos_x.push_back(1.0 + i * 0.5);
            cfg.pos_y.push_back(1.0);
            cfg.pos_z.push_back(1.0);
            cfg.vel_x.push_back(0.1);
            cfg.vel_y.push_back(0.0);
            cfg.vel_z.push_back(0.0);
        }

        auto engine = build_simulation(cfg);
        engine.run();
    }

    EXPECT_TRUE(std::filesystem::exists(ckpt_path));

    {
        SimulationConfig cfg;
        cfg.dt = 0.002;
        cfg.n_steps = 10;
        cfg.box_size = 3.0;
        cfg.use_lj = true;
        cfg.lj = LennardJonesParams{};
        cfg.lj.n_types = 1;
        cfg.lj.sigma = {0.3405};
        cfg.lj.epsilon = {0.997};
        cfg.lj_cutoff = 1.2;

        int n = 4;
        cfg.masses.assign(n, 39.948);
        cfg.charges.assign(n, 0.0);
        for (int i = 0; i < n; ++i) {
            cfg.pos_x.push_back(1.0 + i * 0.5);
            cfg.pos_y.push_back(1.0);
            cfg.pos_z.push_back(1.0);
            cfg.vel_x.push_back(0.1);
            cfg.vel_y.push_back(0.0);
            cfg.vel_z.push_back(0.0);
        }

        auto engine = build_simulation(cfg);
        engine.load_checkpoint(ckpt_path);
        EXPECT_EQ(engine.current_step(), 5);
        engine.run();
        EXPECT_EQ(engine.current_step(), 15);
    }

    std::filesystem::remove(ckpt_path);
}
