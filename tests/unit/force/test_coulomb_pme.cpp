#include "force/coulomb_pme.h"
#include "core/cell.h"
#include <gtest/gtest.h>
#include <cmath>
#include <random>

// ---- B-spline unit tests ----

TEST(CoulombPMETest, BSplineValues) {
    EXPECT_NEAR(CoulombPME::theta4(0.0), 2.0 / 3.0, 1e-12);
    EXPECT_NEAR(CoulombPME::theta4(1.0), 1.0 / 6.0, 1e-12);
    EXPECT_NEAR(CoulombPME::theta4(2.0), 0.0, 1e-12);
    EXPECT_NEAR(CoulombPME::theta4(0.5), (4.0 - 6.0 * 0.25 + 3.0 * 0.125) / 6.0, 1e-12);
    EXPECT_NEAR(CoulombPME::theta4(0.3), CoulombPME::theta4(-0.3), 1e-12);
}

TEST(CoulombPMETest, BSplineDerivative) {
    EXPECT_NEAR(CoulombPME::dtheta4(0.0), 0.0, 1e-12);
    EXPECT_NEAR(CoulombPME::dtheta4(0.5), -CoulombPME::dtheta4(-0.5), 1e-12);
    EXPECT_NEAR(CoulombPME::dtheta4(1.0), -0.5, 1e-12);
    EXPECT_NEAR(CoulombPME::dtheta4(-1.0), 0.5, 1e-12);
}

// ---- Self-energy tests ----

TEST(CoulombPMETest, SelfEnergySingleCharge) {
    double k_e = 138.935456;
    double alpha = 3.5;
    std::vector<double> charges = {1.0};
    double expected = -k_e * alpha / std::sqrt(M_PI) * 1.0;
    EXPECT_NEAR(CoulombPME::calc_self_energy(alpha, charges, k_e), expected, 1e-10);
}

TEST(CoulombPMETest, SelfEnergyMultipleCharges) {
    double k_e = 138.935456;
    double alpha = 3.5;
    std::vector<double> charges = {1.0, -1.0, 0.5, -0.5};
    double sum_q2 = 1.0 + 1.0 + 0.25 + 0.25;
    double expected = -k_e * alpha / std::sqrt(M_PI) * sum_q2;
    EXPECT_NEAR(CoulombPME::calc_self_energy(alpha, charges, k_e), expected, 1e-10);
}

// ---- Force symmetry: Newton's 3rd law within grid discretization ----

TEST(CoulombPMETest, ForceSymmetry) {
    PMEParams params;
    params.cutoff = 2.0;
    params.ewald_coefficient = 3.5;
    params.nx = 32;
    params.ny = 32;
    params.nz = 32;
    params.spline_order = 4;

    CoulombPME pme(params);

    double L = 4.0;
    Cell cell(L, L, L);

    int n_atoms = 4;
    std::vector<double> pos_x = {1.2, 2.8, 1.2, 2.8};
    std::vector<double> pos_y = {1.2, 1.2, 2.8, 2.8};
    std::vector<double> pos_z = {1.2, 1.2, 2.8, 2.8};
    std::vector<double> charges = {1.0, -1.0, 1.0, -1.0};
    std::vector<double> forces_x(n_atoms, 0.0);
    std::vector<double> forces_y(n_atoms, 0.0);
    std::vector<double> forces_z(n_atoms, 0.0);

    pme.set_charges(charges);

    double energy = 0.0;
    pme.compute(pos_x, pos_y, pos_z, cell,
                forces_x, forces_y, forces_z, energy);

    double net_x = 0, net_y = 0, net_z = 0;
    for (int i = 0; i < n_atoms; ++i) {
        net_x += forces_x[i];
        net_y += forces_y[i];
        net_z += forces_z[i];
    }
    EXPECT_NEAR(net_x, 0.0, 0.02);
    EXPECT_NEAR(net_y, 0.0, 0.02);
    EXPECT_NEAR(net_z, 0.0, 0.02);
}

// ---- Energy conservation for charged system (short NVE) ----

TEST(CoulombPMETest, ShortNVEEnergyDrift) {
    PMEParams params;
    params.cutoff = 2.0;
    params.ewald_coefficient = 3.5;
    params.nx = 32;
    params.ny = 32;
    params.nz = 32;
    params.spline_order = 4;

    CoulombPME pme(params);

    double L = 3.0;
    Cell cell(L, L, L);

    int n = 8;
    std::vector<double> pos_x, pos_y, pos_z, charges;
    std::mt19937_64 rng(42);
    std::uniform_real_distribution<double> dist(0.5, L - 0.5);

    for (int i = 0; i < n; ++i) {
        pos_x.push_back(dist(rng));
        pos_y.push_back(dist(rng));
        pos_z.push_back(dist(rng));
        charges.push_back((i % 2 == 0) ? 1.0 : -1.0);
    }
    pme.set_charges(charges);

    std::vector<double> forces_x(n, 0.0), forces_y(n, 0.0), forces_z(n, 0.0);
    std::vector<double> vel_x(n, 0.0), vel_y(n, 0.0), vel_z(n, 0.0);

    double mass = 39.948;
    double dt = 0.001;
    int n_steps = 50;

    for (int i = 0; i < n; ++i) {
        vel_x[i] = 0.1;
    }

    double e_pot = 0.0;
    pme.compute(pos_x, pos_y, pos_z, cell,
                forces_x, forces_y, forces_z, e_pot);

    double e_kin = 0.0;
    for (int i = 0; i < n; ++i) {
        e_kin += 0.5 * mass * (vel_x[i]*vel_x[i] + vel_y[i]*vel_y[i] + vel_z[i]*vel_z[i]);
    }
    double e_tot0 = e_kin + e_pot;
    double max_drift = 0.0;

    for (int step = 0; step < n_steps; ++step) {
        for (int i = 0; i < n; ++i) {
            vel_x[i] += 0.5 * forces_x[i] / mass * dt;
            vel_y[i] += 0.5 * forces_y[i] / mass * dt;
            vel_z[i] += 0.5 * forces_z[i] / mass * dt;
        }
        for (int i = 0; i < n; ++i) {
            pos_x[i] += vel_x[i] * dt;
            pos_y[i] += vel_y[i] * dt;
            pos_z[i] += vel_z[i] * dt;
            if (pos_x[i] < 0) pos_x[i] += L;
            if (pos_x[i] >= L) pos_x[i] -= L;
            if (pos_y[i] < 0) pos_y[i] += L;
            if (pos_y[i] >= L) pos_y[i] -= L;
            if (pos_z[i] < 0) pos_z[i] += L;
            if (pos_z[i] >= L) pos_z[i] -= L;
        }

        std::fill(forces_x.begin(), forces_x.end(), 0.0);
        std::fill(forces_y.begin(), forces_y.end(), 0.0);
        std::fill(forces_z.begin(), forces_z.end(), 0.0);
        e_pot = 0.0;
        pme.compute(pos_x, pos_y, pos_z, cell,
                    forces_x, forces_y, forces_z, e_pot);

        for (int i = 0; i < n; ++i) {
            vel_x[i] += 0.5 * forces_x[i] / mass * dt;
            vel_y[i] += 0.5 * forces_y[i] / mass * dt;
            vel_z[i] += 0.5 * forces_z[i] / mass * dt;
        }

        e_kin = 0.0;
        for (int i = 0; i < n; ++i) {
            e_kin += 0.5 * mass * (vel_x[i]*vel_x[i] + vel_y[i]*vel_y[i] + vel_z[i]*vel_z[i]);
        }
        double e_tot = e_kin + e_pot;
        double drift = std::abs(e_tot - e_tot0) / (std::abs(e_tot0) + 1e-10);
        if (drift > max_drift) max_drift = drift;
    }
    EXPECT_LT(max_drift, 0.01);
}
