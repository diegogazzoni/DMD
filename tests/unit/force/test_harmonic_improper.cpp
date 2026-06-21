#include "force/harmonic_improper.h"
#include "core/cell.h"
#include "core/system_data.h"
#include <gtest/gtest.h>
#include <cmath>

TEST(HarmonicImproperTest, PlanarEquilibrium) {
    ImproperParams params;
    params.n = 1;
    params.i = {0}; params.j = {1}; params.k = {2}; params.l = {3};
    params.k_phi = {100.0};
    params.phi0 = {M_PI};  // planar: 180°

    HarmonicImproper imp(std::move(params));

    SystemData sys(4);
    // Planar cis configuration: all atoms in z=0 plane, angle around j-k bond = 180°
    sys.pos_x[0] = 0.0; sys.pos_y[0] = 1.0; sys.pos_z[0] = 0.0;  // i
    sys.pos_x[1] = 0.0; sys.pos_y[1] = 0.0; sys.pos_z[1] = 0.0;  // j
    sys.pos_x[2] = 1.0; sys.pos_y[2] = 0.0; sys.pos_z[2] = 0.0;  // k
    sys.pos_x[3] = 2.0; sys.pos_y[3] = 0.0; sys.pos_z[3] = 0.0;  // l

    Cell cell(10.0, 10.0, 10.0);
    double energy = 0.0;
    imp.compute(sys.x(), sys.y(), sys.z(), cell, sys.fx(), sys.fy(), sys.fz(), energy);

    EXPECT_NEAR(energy, 0.0, 1e-12);
}

TEST(HarmonicImproperTest, NonPlanarEnergy) {
    ImproperParams params;
    params.n = 1;
    params.i = {0}; params.j = {1}; params.k = {2}; params.l = {3};
    params.k_phi = {100.0};
    params.phi0 = {M_PI};

    HarmonicImproper imp(std::move(params));

    SystemData sys(4);
    sys.pos_x[0] = 0.0; sys.pos_y[0] = 1.0; sys.pos_z[0] = 0.0;
    sys.pos_x[1] = 0.0; sys.pos_y[1] = 0.0; sys.pos_z[1] = 0.0;
    sys.pos_x[2] = 1.0; sys.pos_y[2] = 0.0; sys.pos_z[2] = 0.0;
    sys.pos_x[3] = 2.0; sys.pos_y[3] = 0.0; sys.pos_z[3] = 0.5;  // l out of plane

    Cell cell(10.0, 10.0, 10.0);
    double energy = 0.0;
    imp.compute(sys.x(), sys.y(), sys.z(), cell, sys.fx(), sys.fy(), sys.fz(), energy);

    EXPECT_GT(energy, 0.0);
    EXPECT_LT(energy, 200.0);
}

TEST(HarmonicImproperTest, ForcesSumToZero) {
    ImproperParams params;
    params.n = 1;
    params.i = {0}; params.j = {1}; params.k = {2}; params.l = {3};
    params.k_phi = {100.0};
    params.phi0 = {M_PI};

    HarmonicImproper imp(std::move(params));

    SystemData sys(4);
    sys.pos_x[0] = 0.0; sys.pos_y[0] = 1.0; sys.pos_z[0] = 0.5;
    sys.pos_x[1] = 0.0; sys.pos_y[1] = 0.0; sys.pos_z[1] = 0.0;
    sys.pos_x[2] = 1.0; sys.pos_y[2] = 0.0; sys.pos_z[2] = 0.0;
    sys.pos_x[3] = 1.5; sys.pos_y[3] = 0.5; sys.pos_z[3] = 0.3;

    Cell cell(10.0, 10.0, 10.0);
    double energy = 0.0;
    imp.compute(sys.x(), sys.y(), sys.z(), cell, sys.fx(), sys.fy(), sys.fz(), energy);

    double fx = 0, fy = 0, fz = 0;
    for (int i = 0; i < 4; ++i) {
        fx += sys.forces_x[i];
        fy += sys.forces_y[i];
        fz += sys.forces_z[i];
    }

    EXPECT_NEAR(fx, 0.0, 1e-12);
    EXPECT_NEAR(fy, 0.0, 1e-12);
    EXPECT_NEAR(fz, 0.0, 1e-12);
}

TEST(HarmonicImproperTest, MultipleImpropers) {
    ImproperParams params;
    params.n = 2;
    params.i = {0, 4}; params.j = {1, 5}; params.k = {2, 6}; params.l = {3, 7};
    params.k_phi = {50.0, 200.0};
    params.phi0 = {M_PI, 0.0};

    HarmonicImproper imp(std::move(params));

    SystemData sys(8);
    // First improper: planar, zero energy
    sys.pos_x[0] = 0.0; sys.pos_y[0] = 1.0; sys.pos_z[0] = 0.0;
    sys.pos_x[1] = 0.0; sys.pos_y[1] = 0.0; sys.pos_z[1] = 0.0;
    sys.pos_x[2] = 1.0; sys.pos_y[2] = 0.0; sys.pos_z[2] = 0.0;
    sys.pos_x[3] = 2.0; sys.pos_y[3] = 0.0; sys.pos_z[3] = 0.0;
    // Second improper: phi0=0, at 90° → energy
    sys.pos_x[4] = 0.0; sys.pos_y[4] = 1.0; sys.pos_z[4] = 0.0;
    sys.pos_x[5] = 0.0; sys.pos_y[5] = 0.0; sys.pos_z[5] = 0.0;
    sys.pos_x[6] = 1.0; sys.pos_y[6] = 0.0; sys.pos_z[6] = 0.0;
    sys.pos_x[7] = 2.0; sys.pos_y[7] = 0.0; sys.pos_z[7] = 0.0;

    Cell cell(10.0, 10.0, 10.0);
    double energy = 0.0;
    imp.compute(sys.x(), sys.y(), sys.z(), cell, sys.fx(), sys.fy(), sys.fz(), energy);

    EXPECT_NEAR(energy, 0.0, 1e-10);

    // Move atom 7 up to create chi=90°
    sys.pos_z[7] = 1.0;
    energy = 0.0;
    std::fill(sys.forces_x.begin(), sys.forces_x.end(), 0.0);
    std::fill(sys.forces_y.begin(), sys.forces_y.end(), 0.0);
    std::fill(sys.forces_z.begin(), sys.forces_z.end(), 0.0);
    imp.compute(sys.x(), sys.y(), sys.z(), cell, sys.fx(), sys.fy(), sys.fz(), energy);

    // chi = 90°, chi0 = 0°, k = 200 → E = 0.5 * 200 * (π/2)² ≈ 246.74
    double expected = 0.5 * 200.0 * (M_PI / 2.0) * (M_PI / 2.0);
    EXPECT_NEAR(energy, expected, 1e-10);
}
