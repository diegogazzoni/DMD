#include "barostat/barostat.h"
#include "barostat/berendsen_barostat.h"
#include "barostat/andersen_barostat.h"
#include "core/system_data.h"
#include "core/cell.h"
#include <gtest/gtest.h>
#include <cmath>

TEST(BarostatTest, VirialAndPressure) {
    int n = 10;
    SystemData sys(n);
    Cell cell(3.0, 3.0, 3.0);

    for (int i = 0; i < n; ++i) {
        sys.masses[i] = 1.0;
        sys.vel_x[i] = 0.0;
        sys.vel_y[i] = 0.0;
        sys.vel_z[i] = 0.0;
        sys.pos_x[i] = static_cast<double>(i);
        sys.pos_y[i] = 0.0;
        sys.pos_z[i] = 0.0;
        sys.forces_x[i] = 0.1;
        sys.forces_y[i] = 0.0;
        sys.forces_z[i] = 0.0;
    }

    double w = virial(sys);
    double expected_w = 0.0;
    for (int i = 0; i < n; ++i) {
        expected_w += static_cast<double>(i) * 0.1;
    }
    EXPECT_NEAR(w, expected_w, 1e-12);

    double P = pressure(sys, cell.volume());
    EXPECT_GT(P, 0);
}

TEST(BarostatTest, BerendsenBarostatScales) {
    int n = 100;
    SystemData sys(n);
    Cell cell(3.0, 3.0, 3.0);

    for (int i = 0; i < n; ++i) {
        sys.masses[i] = 39.948;
        sys.vel_x[i] = 0.5;
        sys.vel_y[i] = 0.0;
        sys.vel_z[i] = 0.0;
        sys.pos_x[i] = 1.0 + static_cast<double>(i) * 0.01;
        sys.pos_y[i] = 1.0;
        sys.pos_z[i] = 1.0;
        sys.forces_x[i] = 0.0;
        sys.forces_y[i] = 0.0;
        sys.forces_z[i] = 0.0;
    }

    BerendsenBarostatParams params;
    params.pressure = 1.0; // 1 bar
    params.tau = 0.1; // strong coupling
    BerendsenBarostat baro(params);

    double vol_before = cell.volume();
    baro.apply(sys, cell, 0.002);
    double vol_after = cell.volume();

    EXPECT_NE(vol_before, vol_after);
}

TEST(BarostatTest, AndersenBarostatChangesVolume) {
    int n = 100;
    SystemData sys(n);
    Cell cell(3.0, 3.0, 3.0);

    for (int i = 0; i < n; ++i) {
        sys.masses[i] = 39.948;
        sys.vel_x[i] = 0.5;
        sys.vel_y[i] = 0.0;
        sys.vel_z[i] = 0.0;
        sys.pos_x[i] = 1.0 + static_cast<double>(i) * 0.01;
        sys.pos_y[i] = 1.0;
        sys.pos_z[i] = 1.0;
        sys.forces_x[i] = 0.0;
        sys.forces_y[i] = 0.0;
        sys.forces_z[i] = 0.0;
    }

    AndersenBarostatParams params;
    params.pressure = 1.0;
    params.tau = 0.1;
    AndersenBarostat baro(params, n);

    double vol_before = cell.volume();
    baro.apply(sys, cell, 0.002);
    double vol_after = cell.volume();

    EXPECT_NE(vol_before, vol_after);
}

TEST(BarostatTest, PressureZeroVolume) {
    SystemData sys(1);
    EXPECT_DOUBLE_EQ(pressure(sys, 0.0), 0.0);
}
