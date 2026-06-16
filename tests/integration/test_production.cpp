#include "sim/simulation_config.h"
#include "sim/simulation_engine.h"
#include "sim/json_config.h"
#include "sysbin/dmdin.h"
#include "core/cell.h"
#include <gtest/gtest.h>
#include <hdf5.h>
#include <filesystem>
#include <fstream>
#include <vector>

static constexpr double AR_SIGMA = 0.3405;
static constexpr double AR_EPSILON = 0.997;
static constexpr double AR_MASS = 39.948;

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

static std::string make_config_json() {
    return R"({
        "run": {
            "dt": 0.002,
            "n_steps": 50,
            "init_temperature": 85.0,
            "gen_vel": false,
            "seed": 42
        },
        "output": {
            "trajectory_path": "prod_traj.h5",
            "nstxout": 10,
            "nstvout": 10,
            "nstenergy": 5,
            "energy_path": "prod_energy.h5",
            "checkpoint_interval": 0,
            "checkpoint_path": "prod_ckpt.bin"
        },
        "lj": {
            "cutoff": 1.0
        },
        "electrostatics": {
            "coulomb_type": "cutoff",
            "cutoff": 1.0,
            "pme_order": 4,
            "pme_grid_spacing": 0.12,
            "ewald_coeff": 3.5
        },
        "thermostat": {
            "type": "none",
            "temperature": 85.0,
            "tau": 0.1,
            "frequency": 10.0
        },
        "barostat": {
            "type": "none",
            "pressure": 1.0,
            "tau": 1.0,
            "compressibility": 4.5e-5
        },
        "constraints": {
            "type": "none",
            "tolerance": 1e-6
        }
    })";
}

TEST(ProductionGate, EndToEndDmdinConfigToH5MD) {
    std::string dmdin_path = "prod_test_system.dmdin";
    std::string config_path = "prod_test_config.json";
    std::string traj_path = "prod_traj.h5";

    // Clean up any leftovers
    std::filesystem::remove(dmdin_path);
    std::filesystem::remove(config_path);
    std::filesystem::remove(traj_path);

    // 1. Create and write input system
    {
        int n_cells = 2;
        size_t n_atoms = 4 * n_cells * n_cells * n_cells;
        double box_nm = 1.5;

        SimulationConfig cfg;
        cfg.box_size = box_nm;
        cfg.use_lj = true;
        cfg.lj = LennardJonesParams{};
        cfg.lj.n_types = 1;
        cfg.lj.sigma = {AR_SIGMA};
        cfg.lj.epsilon = {AR_EPSILON};
        cfg.lj_cutoff = 1.0;
        cfg.masses.assign(n_atoms, AR_MASS);
        cfg.charges.assign(n_atoms, 0.0);

        fcc_lattice(cfg.pos_x, cfg.pos_y, cfg.pos_z, box_nm, n_cells);
        cfg.vel_x.assign(n_atoms, 0.0);
        cfg.vel_y.assign(n_atoms, 0.0);
        cfg.vel_z.assign(n_atoms, 0.0);

        write_dmdin(dmdin_path, cfg);
    }

    // 2. Write config JSON
    {
        std::ofstream os(config_path);
        ASSERT_TRUE(os);
        os << make_config_json();
    }

    // 3. Full pipeline: read → config → build → run
    {
        SimulationConfig cfg = read_dmdin(dmdin_path);
        apply_json_config(cfg, config_path);

        SimulationEngine engine = build_simulation(cfg);
        engine.run();

        EXPECT_EQ(engine.current_step(), 50);
    }

    // 4. Verify trajectory H5MD exists and is valid
    ASSERT_TRUE(std::filesystem::exists(traj_path));

    hid_t file = H5Fopen(traj_path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    ASSERT_GE(file, 0);

    // Check version
    hid_t ds = H5Dopen2(file, "/h5md/version", H5P_DEFAULT);
    ASSERT_GE(ds, 0);
    int version[2];
    H5Dread(ds, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, version);
    EXPECT_EQ(version[0], 1);
    EXPECT_EQ(version[1], 1);
    H5Dclose(ds);

    // Check frame count: 50 steps, nstxout=10 => 5 frames
    ds = H5Dopen2(file, "/particles/atoms/position/value", H5P_DEFAULT);
    ASSERT_GE(ds, 0);
    hid_t space = H5Dget_space(ds);
    hsize_t dims[3];
    H5Sget_simple_extent_dims(space, dims, nullptr);
    EXPECT_EQ(dims[0], 5);  // 50/10 = 5 frames
    EXPECT_EQ(dims[2], 3);  // 3 spatial dimensions
    H5Sclose(space);
    H5Dclose(ds);

    H5Fclose(file);

    // Cleanup
    std::filesystem::remove(dmdin_path);
    std::filesystem::remove(config_path);
    std::filesystem::remove(traj_path);
}
