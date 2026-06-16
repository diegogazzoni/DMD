#include "force/lennard_jones.h"
#include "core/cell.h"
#include "core/system_data.h"
#include <gtest/gtest.h>
#include <cmath>

TEST(LennardJonesTest, TwoAtoms) {
    LennardJonesParams params;
    params.n_types = 1;
    params.sigma = {0.3, 0.3, 0.3, 0.3}; // [1][1] row-major = 0.336*0.336=0.112...
    params.epsilon = {0.5, 0.5, 0.5, 0.5}; // 1.0 kJ/mol

    double sig = 0.3;
    double eps = 0.5;
    params.sigma = {sig * sig};
    params.epsilon = {eps};
    params.n_types = 1;

    LennardJones lj(std::move(params), 2.0);

    SystemData sys(2);
    sys.pos_x[0] = 0.0; sys.pos_y[0] = 0.0; sys.pos_z[0] = 0.0;
    sys.pos_x[1] = sig; sys.pos_y[1] = 0.0; sys.pos_z[1] = 0.0;
    sys.atom_types[0] = 0;
    sys.atom_types[1] = 0;

    Cell cell(10.0, 10.0, 10.0);
    lj.rebuild_pair_list(sys.x(), sys.y(), sys.z(), cell);

    double energy = 0.0;
    lj.compute(sys.x(), sys.y(), sys.z(), cell, sys.fx(), sys.fy(), sys.fz(), energy);

    EXPECT_NEAR(energy, 0.0, 1e-12);
}

TEST(LennardJonesTest, TwoAtomsRepulsive) {
    LennardJonesParams params;
    params.n_types = 1;
    params.sigma = {0.2649};
    params.epsilon = {0.5};

    LennardJones lj(std::move(params), 2.0);

    SystemData sys(2);
    sys.pos_x[0] = 0.0; sys.pos_y[0] = 0.0; sys.pos_z[0] = 0.0;
    sys.pos_x[1] = 0.2; sys.pos_y[1] = 0.0; sys.pos_z[1] = 0.0;
    sys.atom_types[0] = 0;
    sys.atom_types[1] = 0;

    Cell cell(10.0, 10.0, 10.0);
    lj.rebuild_pair_list(sys.x(), sys.y(), sys.z(), cell);

    double energy = 0.0;
    lj.compute(sys.x(), sys.y(), sys.z(), cell, sys.fx(), sys.fy(), sys.fz(), energy);

    EXPECT_GT(energy, 0.0);
    EXPECT_GT(sys.forces_x[0], 0.0); // repulsive: atom0 pushed away from atom1
    EXPECT_LT(sys.forces_x[1], 0.0);
}
