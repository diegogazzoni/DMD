#include "core/system_data.h"
#include <gtest/gtest.h>

TEST(SystemDataTest, DefaultConstruction) {
    SystemData sys;
    EXPECT_EQ(sys.n_atoms, 0);
    EXPECT_DOUBLE_EQ(sys.potential_energy, 0.0);
}

TEST(SystemDataTest, SizedConstruction) {
    SystemData sys(100);
    EXPECT_EQ(sys.n_atoms, 100);
    EXPECT_EQ(sys.x().size(), 100);
    EXPECT_EQ(sys.fx().size(), 100);
    EXPECT_EQ(sys.masses.size(), 100);
}

TEST(SystemDataTest, SpanAccess) {
    SystemData sys(10);
    sys.pos_x[3] = 1.5;
    sys.pos_y[3] = 2.5;

    EXPECT_DOUBLE_EQ(sys.x()[3], 1.5);
    EXPECT_DOUBLE_EQ(sys.y()[3], 2.5);
}

TEST(SystemDataTest, ForceSpanMutable) {
    SystemData sys(5);
    sys.fx()[2] = -3.0;
    EXPECT_DOUBLE_EQ(sys.forces_x[2], -3.0);
}
