#include "force/force_engine.h"
#include "force/lennard_jones.h"
#include "core/cell.h"
#include "core/system_data.h"
#include "integrate/integrator.h"
#include <gtest/gtest.h>
#include <cmath>
#include <random>
#include <vector>

static constexpr double AR_SIGMA = 0.3405;   // nm
static constexpr double AR_EPSILON = 0.997;  // kJ/mol
static constexpr double AR_MASS = 39.948;    // amu
static constexpr double KB = 0.008314462618; // kJ/(mol*K)

static void maxwell_boltzmann(
    std::span<double> vx,
    std::span<double> vy,
    std::span<double> vz,
    std::span<const double> masses,
    double temperature)
{
    std::mt19937_64 rng(42);
    std::normal_distribution<double> dist(0.0, 1.0);

    for (size_t i = 0; i < vx.size(); ++i) {
        double sigma_v = std::sqrt(KB * temperature / masses[i]);
        vx[i] = dist(rng) * sigma_v;
        vy[i] = dist(rng) * sigma_v;
        vz[i] = dist(rng) * sigma_v;
    }

    // remove center-of-mass motion
    double com_vx = 0, com_vy = 0, com_vz = 0, total_mass = 0;
    for (size_t i = 0; i < vx.size(); ++i) {
        com_vx += masses[i] * vx[i];
        com_vy += masses[i] * vy[i];
        com_vz += masses[i] * vz[i];
        total_mass += masses[i];
    }
    com_vx /= total_mass;
    com_vy /= total_mass;
    com_vz /= total_mass;

    for (size_t i = 0; i < vx.size(); ++i) {
        vx[i] -= com_vx;
        vy[i] -= com_vy;
        vz[i] -= com_vz;
    }

    // rescale to target temperature
    double ekin = 0;
    for (size_t i = 0; i < vx.size(); ++i) {
        ekin += 0.5 * masses[i] * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
    }
    double dof = 3.0 * vx.size() - 3.0;
    double T_actual = 2.0 * ekin / (dof * KB);
    double scale = std::sqrt(temperature / T_actual);
    for (size_t i = 0; i < vx.size(); ++i) {
        vx[i] *= scale;
        vy[i] *= scale;
        vz[i] *= scale;
    }
}

static void fcc_lattice(
    std::span<double> x,
    std::span<double> y,
    std::span<double> z,
    double box_size,
    int n_cells)
{
    // 4 atoms per FCC unit cell
    double a = box_size / n_cells;
    int idx = 0;
    for (int i = 0; i < n_cells; ++i) {
        for (int j = 0; j < n_cells; ++j) {
            for (int k = 0; k < n_cells; ++k) {
                double base_x = i * a;
                double base_y = j * a;
                double base_z = k * a;

                x[idx] = base_x;          y[idx] = base_y;          z[idx] = base_z;          ++idx;
                x[idx] = base_x + a*0.5;  y[idx] = base_y + a*0.5;  z[idx] = base_z;          ++idx;
                x[idx] = base_x + a*0.5;  y[idx] = base_y;          z[idx] = base_z + a*0.5;  ++idx;
                x[idx] = base_x;          y[idx] = base_y + a*0.5;  z[idx] = base_z + a*0.5;  ++idx;
            }
        }
    }
}

TEST(ArgonNVE, EnergyConservationShort) {
    // 4x4x4 FCC = 256 atoms
    int n_cells = 4;
    size_t n_atoms = 4 * n_cells * n_cells * n_cells;
    double density = 1.4; // g/cm^3 -> ~0.84 g/mL for liquid argon
    // density = mass/volume -> volume = mass/density
    // mass per atom in g: 39.948 / 6.022e23 = 6.634e-23 g
    // for n_atoms: total_mass_g = n_atoms * 6.634e-23 g
    // volume = total_mass_g / density (g/cm^3) = in cm^3
    // box length = volume^(1/3) in cm -> convert to nm: * 1e7
    // box length in nm for density 1.4 g/cm^3:
    double mass_total_g = n_atoms * AR_MASS / 6.02214076e23;
    double volume_cm3 = mass_total_g / density;
    double box_nm = std::pow(volume_cm3, 1.0/3.0) * 1e7;

    SystemData sys(n_atoms);
    sys.masses.assign(n_atoms, AR_MASS);

    fcc_lattice(sys.pos_x, sys.pos_y, sys.pos_z, box_nm, n_cells);

    Cell cell(box_nm, box_nm, box_nm);

    double temperature = 85.0; // K
    maxwell_boltzmann(sys.vel_x, sys.vel_y, sys.vel_z, sys.masses, temperature);

    // set up LJ
    LennardJonesParams lj_params;
    lj_params.n_types = 1;
    lj_params.sigma = {AR_SIGMA};
    lj_params.epsilon = {AR_EPSILON};

    auto lj = std::make_unique<LennardJones>(std::move(lj_params), box_nm * 0.49);
    lj->set_atom_types(sys.atom_types);

    // set up ForceEngine
    ForceEngine engine;
    engine.add_component(std::move(lj));

    // initial forces
    sys.potential_energy = 0.0;
    engine.compute(sys, cell);

    Integrator integrator;
    double dt = 0.002; // ps

    // initial total energy
    double ekin0 = 0;
    for (size_t i = 0; i < n_atoms; ++i) {
        ekin0 += 0.5 * AR_MASS * (
            sys.vel_x[i]*sys.vel_x[i] +
            sys.vel_y[i]*sys.vel_y[i] +
            sys.vel_z[i]*sys.vel_z[i]);
    }
    double etot0 = ekin0 + sys.potential_energy;

    int n_steps = 5000;
    double ekin, etot;
    double max_drift = 0.0;

    for (int step = 0; step < n_steps; ++step) {
        integrator.half_kick(sys, dt);
        integrator.advance(sys, dt);
        sys.forces_x.assign(n_atoms, 0.0);
        sys.forces_y.assign(n_atoms, 0.0);
        sys.forces_z.assign(n_atoms, 0.0);
        sys.potential_energy = 0.0;
        engine.compute(sys, cell);
        integrator.half_kick(sys, dt);

        ekin = 0;
        for (size_t i = 0; i < n_atoms; ++i) {
            ekin += 0.5 * AR_MASS * (
                sys.vel_x[i]*sys.vel_x[i] +
                sys.vel_y[i]*sys.vel_y[i] +
                sys.vel_z[i]*sys.vel_z[i]);
        }
        etot = ekin + sys.potential_energy;
        double drift = std::abs(etot - etot0) / std::abs(etot0);
        if (drift > max_drift) max_drift = drift;
    }

    // Relative energy drift should be small (< 0.5% over 5000 steps for dt=0.002)
    EXPECT_LT(max_drift, 5e-3);
}
