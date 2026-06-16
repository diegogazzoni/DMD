#include "force/coulomb_direct.h"
#include "core/cell.h"
#include "core/system_data.h"
#include <gtest/gtest.h>
#include <cmath>

TEST(CoulombDirectTest, TwoEqualCharges) {
    CoulombDirectParams params;
    params.cutoff = 2.0;
    params.ewald_coefficient = 0.5;

    CoulombDirect coulomb(params);
    SystemData sys(2);
    sys.pos_x[0] = 0.0; sys.pos_y[0] = 0.0; sys.pos_z[0] = 0.0;
    sys.pos_x[1] = 1.0; sys.pos_y[1] = 0.0; sys.pos_z[1] = 0.0;

    std::vector<double> charges = {1.0, 1.0};
    coulomb.set_charges(charges);

    Cell cell(10.0, 10.0, 10.0);
    coulomb.rebuild_pair_list(sys.x(), sys.y(), sys.z(), cell);

    double energy = 0.0;
    coulomb.compute(sys.x(), sys.y(), sys.z(), cell, sys.fx(), sys.fy(), sys.fz(), energy);

    double k_e = 138.935456;
    double expected_energy = k_e * std::erfc(0.5) / 1.0;
    EXPECT_GT(energy, 0.0);
}

TEST(CoulombDirectTest, OppositeChargesAttractive) {
    CoulombDirectParams params;
    params.cutoff = 2.0;
    params.ewald_coefficient = 0.5;

    CoulombDirect coulomb(params);
    SystemData sys(2);
    sys.pos_x[0] = 0.0; sys.pos_y[0] = 0.0; sys.pos_z[0] = 0.0;
    sys.pos_x[1] = 0.8; sys.pos_y[1] = 0.0; sys.pos_z[1] = 0.0;

    std::vector<double> charges = {1.0, -1.0};
    coulomb.set_charges(charges);

    Cell cell(10.0, 10.0, 10.0);
    coulomb.rebuild_pair_list(sys.x(), sys.y(), sys.z(), cell);

    double energy = 0.0;
    coulomb.compute(sys.x(), sys.y(), sys.z(), cell, sys.fx(), sys.fy(), sys.fz(), energy);

    EXPECT_LT(energy, 0.0);
    EXPECT_GT(sys.forces_x[0], 0.0); // repelled toward center
    EXPECT_LT(sys.forces_x[1], 0.0);
}
