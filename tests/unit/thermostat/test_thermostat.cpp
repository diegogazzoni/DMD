#include "thermostat/thermostat.h"
#include "thermostat/andersen.h"
#include "thermostat/berendsen.h"
#include "thermostat/nose_hoover.h"
#include "core/system_data.h"
#include <gtest/gtest.h>
#include <cmath>

TEST(ThermostatTest, KineticEnergyAndTemperature) {
    int n = 100;
    SystemData sys(n);
    for (int i = 0; i < n; ++i) {
        sys.masses[i] = 1.0;
        sys.vel_x[i] = 0.1;
        sys.vel_y[i] = 0.2;
        sys.vel_z[i] = 0.3;
    }

    double ekin = kinetic_energy(sys);
    double expected_ekin = n * 0.5 * 1.0 * (0.01 + 0.04 + 0.09);
    EXPECT_NEAR(ekin, expected_ekin, 1e-12);

    double T = temperature(ekin, n);
    EXPECT_GT(T, 0);
    EXPECT_EQ(temperature(sys), T);
}

TEST(ThermostatTest, AndersenModifiesVelocities) {
    int n = 1000;
    SystemData sys(n);
    for (int i = 0; i < n; ++i) {
        sys.masses[i] = 39.948;
        sys.vel_x[i] = 0.1;
        sys.vel_y[i] = 0.0;
        sys.vel_z[i] = 0.0;
    }

    AndersenParams params;
    params.temperature = 300.0;
    params.frequency = 500.0; // high freq to ensure many collisions
    AndersenThermostat thermo(params);

    double ekin_before = kinetic_energy(sys);
    thermo.apply(sys, 0.002);
    double ekin_after = kinetic_energy(sys);

    // velocities should have been modified
    EXPECT_NE(ekin_before, ekin_after);
    EXPECT_GT(ekin_after, 0);
}

TEST(ThermostatTest, BerendsenRescalesToTargetT) {
    int n = 1000;
    SystemData sys(n);
    for (int i = 0; i < n; ++i) {
        sys.masses[i] = 39.948;
        sys.vel_x[i] = 0.5;
        sys.vel_y[i] = 0.0;
        sys.vel_z[i] = 0.0;
    }

    BerendsenParams params;
    params.temperature = 100.0; // well below current T
    params.tau = 0.01; // very strong coupling
    BerendsenThermostat thermo(params);

    double T_before = temperature(sys);
    EXPECT_GT(T_before, params.temperature);

    // multiple steps to converge
    for (int i = 0; i < 100; ++i) {
        thermo.apply(sys, 0.002);
    }

    double T_after = temperature(sys);
    EXPECT_NEAR(T_after, params.temperature, params.temperature * 0.05);
}

TEST(ThermostatTest, NoseHooverModifiesVelocities) {
    int n = 100;
    SystemData sys(n);
    for (int i = 0; i < n; ++i) {
        sys.masses[i] = 1.0;
        sys.vel_x[i] = 0.5;
        sys.vel_y[i] = 0.0;
        sys.vel_z[i] = 0.0;
    }

    NoseHooverParams params;
    params.temperature = 100.0;
    params.tau = 0.1;
    NoseHooverThermostat thermo(params, n);

    double ekin_before = kinetic_energy(sys);
    thermo.apply(sys, 0.002);
    double ekin_after = kinetic_energy(sys);

    // zeta should have changed
    EXPECT_NE(thermo.zeta(), 0.0);
    // velocities should be different
    EXPECT_NE(ekin_before, ekin_after);
}

TEST(ThermostatTest, BerendsenZeroTauNoOp) {
    int n = 10;
    SystemData sys(n);
    for (int i = 0; i < n; ++i) {
        sys.masses[i] = 1.0;
        sys.vel_x[i] = 1.0;
        sys.vel_y[i] = 0.0;
        sys.vel_z[i] = 0.0;
    }

    BerendsenParams params;
    params.tau = 0.0;
    BerendsenThermostat thermo(params);

    double ekin_before = kinetic_energy(sys);
    thermo.apply(sys, 0.002);
    double ekin_after = kinetic_energy(sys);

    EXPECT_DOUBLE_EQ(ekin_before, ekin_after);
}
