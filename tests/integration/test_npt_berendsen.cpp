#include "sim/simulation_engine.h"
#include "force/force_engine.h"
#include "force/lennard_jones.h"
#include "thermostat/berendsen.h"
#include "thermostat/thermostat.h"
#include "barostat/berendsen_barostat.h"
#include "barostat/barostat.h"
#include "core/constants.h"
#include <gtest/gtest.h>
#include <cmath>
#include <vector>
#include <random>

static constexpr double AR_SIGMA = 0.3405;
static constexpr double AR_EPSILON = 0.997;
static constexpr double AR_MASS = 39.948;

static void fcc_lattice(std::vector<double>& x, std::vector<double>& y, std::vector<double>& z,
                        double box, int n) {
    double a = box / n;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                double bx = i * a, by = j * a, bz = k * a;
                x.push_back(bx);          y.push_back(by);          z.push_back(bz);
                x.push_back(bx + a*0.5);  y.push_back(by + a*0.5);  z.push_back(bz);
                x.push_back(bx + a*0.5);  y.push_back(by);          z.push_back(bz + a*0.5);
                x.push_back(bx);          y.push_back(by + a*0.5);  z.push_back(bz + a*0.5);
            }
        }
    }
}

TEST(NPTGate, BerendsenNPTArgon) {
    int n_cells = 3;
    size_t n_atoms = 4 * n_cells * n_cells * n_cells; // 108

    double target_T = 120.0;
    double target_P = 100.0; // bar

    // initial density ~1.0 g/cm^3
    double initial_density = 1.0;
    double mass_total_g = n_atoms * AR_MASS / 6.02214076e23;
    double box_nm = std::pow(mass_total_g / initial_density, 1.0 / 3.0) * 1e7;

    // positions
    std::vector<double> px, py, pz;
    fcc_lattice(px, py, pz, box_nm, n_cells);

    // Maxwell-Boltzmann velocities at target_T
    std::mt19937_64 rng(42);
    std::normal_distribution<double> ndist(0.0, 1.0);
    double sigma_v = std::sqrt(dmd_constants::KB * target_T / AR_MASS);
    std::vector<double> vx(n_atoms), vy(n_atoms), vz(n_atoms);
    for (size_t i = 0; i < n_atoms; ++i) {
        vx[i] = ndist(rng) * sigma_v;
        vy[i] = ndist(rng) * sigma_v;
        vz[i] = ndist(rng) * sigma_v;
    }
    double com_vx = 0, com_vy = 0, com_vz = 0;
    for (size_t i = 0; i < n_atoms; ++i) {
        com_vx += vx[i]; com_vy += vy[i]; com_vz += vz[i];
    }
    com_vx /= n_atoms; com_vy /= n_atoms; com_vz /= n_atoms;
    for (size_t i = 0; i < n_atoms; ++i) {
        vx[i] -= com_vx; vy[i] -= com_vy; vz[i] -= com_vz;
    }

    // Force engine: LJ only
    auto fe = std::make_unique<ForceEngine>();
    LennardJonesParams lj_params;
    lj_params.n_types = 1;
    lj_params.sigma = {AR_SIGMA};
    lj_params.epsilon = {AR_EPSILON};
    auto lj = std::make_unique<LennardJones>(lj_params, box_nm * 0.49);
    auto* lj_ptr = lj.get();
    fe->add_component(std::move(lj));

    // Cell
    Cell cell(box_nm, box_nm, box_nm);

    // Engine config
    SimulationEngine::Config cfg;
    cfg.dt = 0.002;
    cfg.n_steps = 2000;
    cfg.checkpoint_interval = 0;

    // Thermostat + Barostat
    BerendsenParams t_params;
    t_params.temperature = target_T;
    t_params.tau = 0.1;
    auto thermo = std::make_unique<BerendsenThermostat>(t_params);

    BerendsenBarostatParams b_params;
    b_params.pressure = target_P;
    b_params.tau = 1.0;
    b_params.compressibility = 4.5e-5; // bar^-1 (approx for liquids)
    auto baro = std::make_unique<BerendsenBarostat>(b_params);

    // Build engine
    SimulationEngine engine(std::move(fe), cell, n_atoms, cfg, std::move(thermo), std::move(baro));
    auto& sys = engine.system();

    // Initialize system
    sys.masses.assign(n_atoms, AR_MASS);
    sys.charges.assign(n_atoms, 0.0);
    sys.pos_x = px; sys.pos_y = py; sys.pos_z = pz;
    sys.vel_x = vx; sys.vel_y = vy; sys.vel_z = vz;
    std::vector<int> types(n_atoms, 0);
    sys.atom_types = types;
    lj_ptr->set_atom_types(sys.atom_types);

    // Run
    engine.run();

    // ---- Checks ----
    double ekin = kinetic_energy(sys);
    double T_final = temperature(ekin, n_atoms);
    double vol = cell.volume();
    double density = mass_total_g / (vol * 1e-21); // g/cm^3

    EXPECT_NEAR(T_final, target_T, target_T * 0.10)
        << "Temperature should be within 10% of target";

    EXPECT_GT(density, 0.9)
        << "Density should remain positive and reasonable for liquid Ar";
    EXPECT_LT(density, 2.0)
        << "Density should not be excessively high";

    // At 120K and 100 bar, liquid Ar density is around 1.1-1.3 g/cm^3
    // Allow wider bounds for this test
    EXPECT_NEAR(density, 1.2, 0.5)
        << "Density should converge toward liquid Ar density";
}
