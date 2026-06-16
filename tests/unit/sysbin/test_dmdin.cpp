#include "sysbin/dmdin.h"
#include "sim/simulation_config.h"
#include "sim/simulation_engine.h"
#include "force/force_engine.h"
#include "core/cell.h"
#include <gtest/gtest.h>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <cstdint>

TEST(DMDinTest, WriteAndReadBack) {
    std::string path = "test_system.dmdin";

    // create a config
    SimulationConfig cfg;
    cfg.dt = 0.002;
    cfg.n_steps = 100;
    cfg.box_size = 3.0;
    cfg.lj_cutoff = 1.2;
    cfg.lj.n_types = 1;
    cfg.lj.sigma = {0.3405};
    cfg.lj.epsilon = {0.997};
    cfg.use_lj = true;
    cfg.checkpoint_path = "test.cpt";

    int n = 8;
    for (int i = 0; i < n; ++i) {
        cfg.pos_x.push_back(0.5 + i * 0.3);
        cfg.pos_y.push_back(0.5 + i * 0.3);
        cfg.pos_z.push_back(0.5 + i * 0.3);
        cfg.vel_x.push_back(0.01 * i);
        cfg.vel_y.push_back(0.0);
        cfg.vel_z.push_back(0.0);
        cfg.masses.push_back(39.948);
        cfg.charges.push_back(0.0);
        cfg.atom_types.push_back(0);
    }

    // write
    write_dmdin(path, cfg);

    // verify file exists and has reasonable size
    EXPECT_TRUE(std::filesystem::exists(path));
    EXPECT_GT(std::filesystem::file_size(path), 200);

    // read back
    SimulationConfig cfg2 = read_dmdin(path);

    // compare — dmdin carries system data only (not runtime params)
    EXPECT_DOUBLE_EQ(cfg2.box_size, 3.0);
    EXPECT_TRUE(cfg2.use_lj);

    EXPECT_EQ(cfg2.pos_x.size(), 8);
    EXPECT_EQ(cfg2.vel_x.size(), 8);
    EXPECT_EQ(cfg2.masses.size(), 8);
    EXPECT_EQ(cfg2.charges.size(), 8);

    for (int i = 0; i < n; ++i) {
        EXPECT_NEAR(cfg2.pos_x[i], 0.5 + i * 0.3, 1e-6);
        EXPECT_NEAR(cfg2.pos_y[i], 0.5 + i * 0.3, 1e-6);
        EXPECT_NEAR(cfg2.pos_z[i], 0.5 + i * 0.3, 1e-6);
    }

    for (int i = 0; i < n; ++i) {
        EXPECT_NEAR(cfg2.vel_x[i], 0.01 * i, 1e-6);
    }

    EXPECT_DOUBLE_EQ(cfg2.masses[0], 39.948);
    EXPECT_DOUBLE_EQ(cfg2.charges[0], 0.0);

    std::filesystem::remove(path);
}

TEST(DMDinTest, RoundTripWithBonds) {
    std::string path = "test_bonded.dmdin";

    SimulationConfig cfg;
    cfg.box_size = 2.0;
    cfg.lj.n_types = 1;
    cfg.lj.sigma = {0.34};
    cfg.lj.epsilon = {1.0};
    cfg.use_lj = true;

    int n = 4;
    for (int i = 0; i < n; ++i) {
        cfg.pos_x.push_back(static_cast<double>(i));
        cfg.pos_y.push_back(0.0);
        cfg.pos_z.push_back(0.0);
        cfg.masses.push_back(12.0);
        cfg.charges.push_back(0.0);
        cfg.atom_types.push_back(0);
    }

    HarmonicBondParams bp;
    bp.n = 3;
    bp.i = {0, 1, 2};
    bp.j = {1, 2, 3};
    bp.k = {300.0, 300.0, 300.0};
    bp.r0 = {0.15, 0.15, 0.15};
    cfg.bonds.push_back(bp);

    HarmonicAngleParams ap;
    ap.n = 1;
    ap.i = {0}; ap.j = {1}; ap.k = {2};
    ap.k_theta = {100.0};
    ap.theta0 = {1.8};
    cfg.angles.push_back(ap);

    write_dmdin(path, cfg);
    SimulationConfig cfg2 = read_dmdin(path);

    EXPECT_EQ(cfg2.bonds.size(), 1);
    EXPECT_EQ(cfg2.bonds[0].n, 3);
    EXPECT_EQ(cfg2.bonds[0].i[0], 0);
    EXPECT_EQ(cfg2.bonds[0].j[1], 2);
    EXPECT_DOUBLE_EQ(cfg2.bonds[0].k[0], 300.0);
    EXPECT_DOUBLE_EQ(cfg2.bonds[0].r0[2], 0.15);

    EXPECT_EQ(cfg2.angles.size(), 1);
    EXPECT_EQ(cfg2.angles[0].n, 1);

    std::filesystem::remove(path);
}

TEST(DMDinTest, ReadSyntheticFileDirectly) {
    // manually construct a minimal dmdin file and read it
    std::string path = "test_synthetic.dmdin";
    std::ofstream os(path, std::ios::binary);

    struct Header {
        uint32_t magic = 0x444D444E;
        uint32_t version = 1;
        uint64_t n_atoms = 2;
        uint32_t cell_type = 0;
        int32_t n_types = 1;
        uint32_t pos_flags = 0;
        uint32_t reserved = 0;
        uint64_t offset_ff = 120 + 2 * 8; // header + masses(1) + charges(1)
        uint64_t offset_pos = 0;
        double box[9] = {2.0,0,0, 0,2.0,0, 0,0,2.0};
    } h;
    h.offset_pos = h.offset_ff + 4 + 4 + 4 + 4 + 8 + 8;
    os.write(reinterpret_cast<const char*>(&h), sizeof(h));

    double mass = 1.0, charge = 0.0;
    os.write(reinterpret_cast<const char*>(&mass), 8);
    os.write(reinterpret_cast<const char*>(&charge), 8);

    // FF section: n_bonds=0, n_angles=0, n_dihedrals=0, n_lj_types=1
    uint32_t zero = 0;
    for (int i = 0; i < 3; ++i) os.write(reinterpret_cast<const char*>(&zero), 4);
    uint32_t lj_n = 1;
    os.write(reinterpret_cast<const char*>(&lj_n), 4);
    double sig = 0.34, eps = 1.0;
    os.write(reinterpret_cast<const char*>(&sig), 8);
    os.write(reinterpret_cast<const char*>(&eps), 8);

    // Positions
    float px[2] = {1.0f, 2.0f}, py[2] = {0,0}, pz[2] = {0,0};
    os.write(reinterpret_cast<const char*>(px), 8);
    os.write(reinterpret_cast<const char*>(py), 8);
    os.write(reinterpret_cast<const char*>(pz), 8);
    int32_t types[2] = {0, 0};
    os.write(reinterpret_cast<const char*>(types), 8);
    os.close();

    auto cfg = read_dmdin(path);
    EXPECT_EQ(cfg.pos_x.size(), 2);
    EXPECT_NEAR(cfg.pos_x[0], 1.0, 1e-6);
    EXPECT_NEAR(cfg.pos_x[1], 2.0, 1e-6);
    EXPECT_NEAR(cfg.masses[0], 1.0, 1e-12);
    EXPECT_EQ(cfg.vel_x.size(), 0); // no velocities

    std::filesystem::remove(path);
}

TEST(DMDinTest, BadMagic) {
    std::string path = "test_bad.dmdin";
    std::ofstream os(path, std::ios::binary);
    uint32_t bad = 0xDEADBEEF;
    os.write(reinterpret_cast<const char*>(&bad), 4);
    os.close();

    EXPECT_THROW(read_dmdin(path), std::runtime_error);
    std::filesystem::remove(path);
}
