#include "force/harmonic_angle.h"
#include "core/cell.h"
#include "core/system_data.h"
#include <gtest/gtest.h>
#include <cmath>

TEST(HarmonicAngleTest, RightAngleEquilibrium) {
    HarmonicAngleParams params;
    params.n = 1;
    params.i = {0}; params.j = {1}; params.k = {2};
    params.k_theta = {100.0};
    params.theta0 = {M_PI / 2.0};

    HarmonicAngle angle(std::move(params));

    SystemData sys(3);
    sys.pos_x[0] = 0.0; sys.pos_y[0] = 0.0;
    sys.pos_x[1] = 0.0; sys.pos_y[1] = 1.0;
    sys.pos_x[2] = 1.0; sys.pos_y[2] = 1.0;

    Cell cell(10.0, 10.0, 10.0);
    double energy = 0.0;
    angle.compute(sys.x(), sys.y(), sys.z(), cell, sys.fx(), sys.fy(), sys.fz(), energy);

    EXPECT_NEAR(energy, 0.0, 1e-12);
}

TEST(HarmonicAngleTest, BentAngle) {
    HarmonicAngleParams params;
    params.n = 1;
    params.i = {0}; params.j = {1}; params.k = {2};
    params.k_theta = {100.0};
    params.theta0 = {M_PI / 2.0};

    HarmonicAngle angle(std::move(params));

    SystemData sys(3);
    sys.pos_x[0] = 0.0; sys.pos_y[0] = 0.0;
    sys.pos_x[1] = 0.0; sys.pos_y[1] = 1.0;
    sys.pos_x[2] = 1.0; sys.pos_y[2] = 2.0;

    Cell cell(10.0, 10.0, 10.0);
    double energy = 0.0;
    angle.compute(sys.x(), sys.y(), sys.z(), cell, sys.fx(), sys.fy(), sys.fz(), energy);

    EXPECT_GT(energy, 0.0);
}

TEST(HarmonicAngleTest, ForcesSumToZero) {
    HarmonicAngleParams params;
    params.n = 1;
    params.i = {0}; params.j = {1}; params.k = {2};
    params.k_theta = {100.0};
    params.theta0 = {M_PI / 3.0};

    HarmonicAngle angle(std::move(params));

    SystemData sys(3);
    sys.pos_x[0] = 0.0; sys.pos_y[0] = 0.0; sys.pos_z[0] = 0.0;
    sys.pos_x[1] = 0.0; sys.pos_y[1] = 1.0; sys.pos_z[1] = 0.0;
    sys.pos_x[2] = 1.0; sys.pos_y[2] = 1.0; sys.pos_z[2] = 1.0;

    Cell cell(10.0, 10.0, 10.0);
    double energy = 0.0;
    angle.compute(sys.x(), sys.y(), sys.z(), cell, sys.fx(), sys.fy(), sys.fz(), energy);

    double fx = sys.forces_x[0] + sys.forces_x[1] + sys.forces_x[2];
    double fy = sys.forces_y[0] + sys.forces_y[1] + sys.forces_y[2];
    double fz = sys.forces_z[0] + sys.forces_z[1] + sys.forces_z[2];

    EXPECT_NEAR(fx, 0.0, 1e-12);
    EXPECT_NEAR(fy, 0.0, 1e-12);
    EXPECT_NEAR(fz, 0.0, 1e-12);
}
