#include "core/cell.h"
#include <gtest/gtest.h>
#include <cmath>

TEST(CellTest, OrthorhombicVolume) {
    Cell cell(2.0, 3.0, 4.0);
    EXPECT_DOUBLE_EQ(cell.volume(), 24.0);
}

TEST(CellTest, MinimumImageOrthorhombic) {
    Cell cell(2.0, 2.0, 2.0);
    double dx = 1.8, dy = 0.3, dz = -0.4;
    double ix, iy, iz;
    cell.minimum_image(dx, dy, dz, ix, iy, iz);
    EXPECT_NEAR(ix, -0.2, 1e-12);
    EXPECT_NEAR(iy, 0.3, 1e-12);
    EXPECT_NEAR(iz, -0.4, 1e-12);
}

TEST(CellTest, Wrap) {
    Cell cell(5.0, 5.0, 5.0);
    double x = 7.3, y = -2.1, z = 5.0;
    cell.wrap(x, y, z);
    EXPECT_NEAR(x, 2.3, 1e-12);
    EXPECT_NEAR(y, 2.9, 1e-12);
    EXPECT_NEAR(z, 0.0, 1e-12);
}

TEST(CellTest, Scale) {
    Cell cell(2.0, 3.0, 4.0);
    cell.scale(1.5, 2.0, 0.5);
    EXPECT_DOUBLE_EQ(cell.matrix[0], 3.0);
    EXPECT_DOUBLE_EQ(cell.matrix[4], 6.0);
    EXPECT_DOUBLE_EQ(cell.matrix[8], 2.0);
}

TEST(CellTest, ZeroVolume) {
    Cell cell(0.0, 5.0, 5.0);
    EXPECT_DOUBLE_EQ(cell.volume(), 0.0);
}
