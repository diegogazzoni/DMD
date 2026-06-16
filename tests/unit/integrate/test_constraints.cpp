#include "integrate/constraints.h"
#include "core/system_data.h"
#include <gtest/gtest.h>
#include <cmath>

TEST(ConstraintsTest, SimpleBondConstraint) {
    ConstraintsParams params;
    params.n = 1;
    params.i = {0}; params.j = {1};
    params.distance_target = {1.0};
    params.algorithm = ConstraintAlgorithm::shake;

    Constraints constraints(std::move(params));

    SystemData sys(2);
    sys.masses[0] = 1.0;
    sys.masses[1] = 1.0;
    sys.pos_x[0] = 0.0; sys.pos_y[0] = 0.0; sys.pos_z[0] = 0.0;
    sys.pos_x[1] = 1.1; sys.pos_y[1] = 0.1; sys.pos_z[1] = 0.05;

    double dt = 0.001;
    constraints.apply(sys, dt);

    double dx = sys.pos_x[1] - sys.pos_x[0];
    double dy = sys.pos_y[1] - sys.pos_y[0];
    double dz = sys.pos_z[1] - sys.pos_z[0];
    double r = std::sqrt(dx*dx + dy*dy + dz*dz);

    EXPECT_NEAR(r, 1.0, 1e-6);
}

TEST(ConstraintsTest, MultipleConstraints) {
    ConstraintsParams params;
    params.n = 2;
    params.i = {0, 1}; params.j = {1, 2};
    params.distance_target = {1.0, 1.0};
    params.algorithm = ConstraintAlgorithm::shake;

    Constraints constraints(std::move(params));

    SystemData sys(3);
    sys.masses[0] = 1.0; sys.masses[1] = 1.0; sys.masses[2] = 1.0;
    sys.pos_x[0] = 0.0;  sys.pos_y[0] = 0.0;  sys.pos_z[0] = 0.0;
    sys.pos_x[1] = 1.05; sys.pos_y[1] = 0.1;  sys.pos_z[1] = 0.0;
    sys.pos_x[2] = 2.1;  sys.pos_y[2] = -0.1; sys.pos_z[2] = 0.0;

    double dt = 0.001;
    constraints.apply(sys, dt);

    double dx01 = sys.pos_x[1] - sys.pos_x[0];
    double dy01 = sys.pos_y[1] - sys.pos_y[0];
    double dz01 = sys.pos_z[1] - sys.pos_z[0];
    double r01 = std::sqrt(dx01*dx01 + dy01*dy01 + dz01*dz01);

    double dx12 = sys.pos_x[2] - sys.pos_x[1];
    double dy12 = sys.pos_y[2] - sys.pos_y[1];
    double dz12 = sys.pos_z[2] - sys.pos_z[1];
    double r12 = std::sqrt(dx12*dx12 + dy12*dy12 + dz12*dz12);

    EXPECT_NEAR(r01, 1.0, 1e-6);
    EXPECT_NEAR(r12, 1.0, 1e-6);
}
