#include "force/periodic_dihedral.h"
#include "core/cell.h"
#include "core/system_data.h"
#include <gtest/gtest.h>
#include <cmath>

TEST(PeriodicDihedralTest, TransPlanarEnergyZero) {
    PeriodicDihedralParams params;
    params.n = 1;
    params.i = {0}; params.j = {1}; params.k = {2}; params.l = {3};
    params.periodicity = {1};
    params.k_phi = {10.0};
    params.phi0 = {M_PI};

    PeriodicDihedral dihedral(std::move(params));

    SystemData sys(4);
    // trans planar: phi = pi
    sys.pos_x[0] = 0.0; sys.pos_y[0] = 0.0; sys.pos_z[0] = 0.0;
    sys.pos_x[1] = 1.0; sys.pos_y[1] = 0.0; sys.pos_z[1] = 0.0;
    sys.pos_x[2] = 2.0; sys.pos_y[2] = 0.0; sys.pos_z[2] = 0.0;
    sys.pos_x[3] = 3.0; sys.pos_y[3] = 0.0; sys.pos_z[3] = 1.0;

    Cell cell(20.0, 20.0, 20.0);
    double energy = 0.0;
    dihedral.compute(sys.x(), sys.y(), sys.z(), cell, sys.fx(), sys.fy(), sys.fz(), energy);

    EXPECT_NEAR(energy, 0.0, 1e-10);
}

TEST(PeriodicDihedralTest, CisPlanar) {
    PeriodicDihedralParams params;
    params.n = 1;
    params.i = {0}; params.j = {1}; params.k = {2}; params.l = {3};
    params.periodicity = {1};
    params.k_phi = {10.0};
    params.phi0 = {0.0};

    PeriodicDihedral dihedral(std::move(params));

    SystemData sys(4);
    // cis planar: phi = 0
    sys.pos_x[0] = 0.0; sys.pos_y[0] = 0.0; sys.pos_z[0] = 0.0;
    sys.pos_x[1] = 1.0; sys.pos_y[1] = 0.0; sys.pos_z[1] = 0.0;
    sys.pos_x[2] = 2.0; sys.pos_y[2] = 0.0; sys.pos_z[2] = 0.0;
    sys.pos_x[3] = 3.0; sys.pos_y[3] = 0.0; sys.pos_z[3] = 0.0;

    Cell cell(20.0, 20.0, 20.0);
    double energy = 0.0;
    dihedral.compute(sys.x(), sys.y(), sys.z(), cell, sys.fx(), sys.fy(), sys.fz(), energy);

    EXPECT_NEAR(energy, 0.0, 1e-10);
}

TEST(PeriodicDihedralTest, ForcesSumToZero) {
    PeriodicDihedralParams params;
    params.n = 1;
    params.i = {0}; params.j = {1}; params.k = {2}; params.l = {3};
    params.periodicity = {3};
    params.k_phi = {5.0};
    params.phi0 = {0.0};

    PeriodicDihedral dihedral(std::move(params));

    SystemData sys(4);
    sys.pos_x[0] = 0.0; sys.pos_y[0] = 0.0; sys.pos_z[0] = 0.0;
    sys.pos_x[1] = 1.0; sys.pos_y[1] = 0.0; sys.pos_z[1] = 0.0;
    sys.pos_x[2] = 1.2; sys.pos_y[2] = 0.5; sys.pos_z[2] = 0.0;
    sys.pos_x[3] = 1.5; sys.pos_y[3] = 0.8; sys.pos_z[3] = 0.5;

    Cell cell(20.0, 20.0, 20.0);
    double energy = 0.0;
    dihedral.compute(sys.x(), sys.y(), sys.z(), cell, sys.fx(), sys.fy(), sys.fz(), energy);

    double fx = sys.forces_x[0] + sys.forces_x[1] + sys.forces_x[2] + sys.forces_x[3];
    double fy = sys.forces_y[0] + sys.forces_y[1] + sys.forces_y[2] + sys.forces_y[3];
    double fz = sys.forces_z[0] + sys.forces_z[1] + sys.forces_z[2] + sys.forces_z[3];

    EXPECT_NEAR(fx, 0.0, 1e-12);
    EXPECT_NEAR(fy, 0.0, 1e-12);
    EXPECT_NEAR(fz, 0.0, 1e-12);
}

TEST(PeriodicDihedralTest, EnergyMatchesAnalytical) {
    PeriodicDihedralParams params;
    params.n = 1;
    params.i = {0}; params.j = {1}; params.k = {2}; params.l = {3};
    params.periodicity = {2};
    params.k_phi = {3.0};
    params.phi0 = {M_PI / 3.0};

    PeriodicDihedral dihedral(std::move(params));

    SystemData sys(4);
    // proper dihedral with non-colinear atoms, phi = pi/2
    // i=(0,0,0) j=(1,0,0) k=(1,1,0) l=(1,1,1)
    sys.pos_x[0] = 0.0; sys.pos_y[0] = 0.0; sys.pos_z[0] = 0.0;
    sys.pos_x[1] = 1.0; sys.pos_y[1] = 0.0; sys.pos_z[1] = 0.0;
    sys.pos_x[2] = 1.0; sys.pos_y[2] = 1.0; sys.pos_z[2] = 0.0;
    sys.pos_x[3] = 1.0; sys.pos_y[3] = 1.0; sys.pos_z[3] = 1.0;

    Cell cell(20.0, 20.0, 20.0);
    double energy = 0.0;
    dihedral.compute(sys.x(), sys.y(), sys.z(), cell, sys.fx(), sys.fy(), sys.fz(), energy);

    double expected = 3.0 * (1.0 + std::cos(2.0 * M_PI / 2.0 - M_PI / 3.0));
    EXPECT_NEAR(energy, expected, 1e-10);
}
