#include "force/harmonic_bond.h"
#include "core/cell.h"
#include "core/system_data.h"
#include <gtest/gtest.h>
#include <cmath>

TEST(HarmonicBondTest, TwoAtomsEquilibrium) {
    HarmonicBondParams params;
    params.n = 1;
    params.i = {0};
    params.j = {1};
    params.k = {100.0};
    params.r0 = {1.0};

    HarmonicBond bond(std::move(params));

    SystemData sys(2);
    sys.pos_x[0] = 0.0; sys.pos_y[0] = 0.0; sys.pos_z[0] = 0.0;
    sys.pos_x[1] = 1.0; sys.pos_y[1] = 0.0; sys.pos_z[1] = 0.0;

    Cell cell(10.0, 10.0, 10.0);
    double energy = 0.0;
    bond.compute(sys.x(), sys.y(), sys.z(), cell, sys.fx(), sys.fy(), sys.fz(), energy);

    EXPECT_NEAR(energy, 0.0, 1e-12);
    EXPECT_NEAR(sys.forces_x[0], 0.0, 1e-12);
    EXPECT_NEAR(sys.forces_x[1], 0.0, 1e-12);
}

TEST(HarmonicBondTest, TwoAtomsStretched) {
    HarmonicBondParams params;
    params.n = 1;
    params.i = {0};
    params.j = {1};
    params.k = {100.0};
    params.r0 = {1.0};

    HarmonicBond bond(std::move(params));

    SystemData sys(2);
    sys.pos_x[0] = 0.0; sys.pos_y[0] = 0.0; sys.pos_z[0] = 0.0;
    sys.pos_x[1] = 1.5; sys.pos_y[1] = 0.0; sys.pos_z[1] = 0.0;

    Cell cell(10.0, 10.0, 10.0);
    double energy = 0.0;
    bond.compute(sys.x(), sys.y(), sys.z(), cell, sys.fx(), sys.fy(), sys.fz(), energy);

    double expected_energy = 0.5 * 100.0 * (0.5 * 0.5);
    EXPECT_NEAR(energy, expected_energy, 1e-12);

    double expected_f = 100.0 * 0.5;
    EXPECT_NEAR(sys.forces_x[0], expected_f, 1e-12);
    EXPECT_NEAR(sys.forces_x[1], -expected_f, 1e-12);
}

TEST(HarmonicBondTest, ThreeBonds) {
    HarmonicBondParams params;
    params.n = 2;
    params.i = {0, 1};
    params.j = {1, 2};
    params.k = {100.0, 200.0};
    params.r0 = {1.0, 1.5};

    HarmonicBond bond(std::move(params));

    SystemData sys(3);
    sys.pos_x[0] = 0.0; sys.pos_y[0] = 0.0;
    sys.pos_x[1] = 1.0; sys.pos_y[1] = 0.0;
    sys.pos_x[2] = 2.8; sys.pos_y[2] = 0.0;

    Cell cell(20.0, 20.0, 20.0);
    double energy = 0.0;
    bond.compute(sys.x(), sys.y(), sys.z(), cell, sys.fx(), sys.fy(), sys.fz(), energy);

    EXPECT_NEAR(energy, 9.0, 1e-10);
}
