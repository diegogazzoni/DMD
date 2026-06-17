#pragma once

#include "core/cell.h"
#include "force/lennard_jones.h"
#include "force/harmonic_bond.h"
#include "force/harmonic_angle.h"
#include "force/periodic_dihedral.h"
#include <memory>
#include <string>
#include <vector>

struct ThermostatConfig {
    std::string type{"none"};
    double temperature{300.0};
    double tau{0.1};
    double frequency{10.0};
};

struct BarostatConfig {
    std::string type{"none"};
    double pressure{1.0};
    double tau{1.0};
    double compressibility{4.5e-5};
};

struct SimulationConfig {
    // Run control
    double dt{0.002};
    int n_steps{1000};
    double init_temperature{300.0};
    bool gen_vel{false};
    int seed{42};

    // Box
    double box_size{3.0};

    // System
    std::vector<double> masses;
    std::vector<double> charges;
    std::vector<double> pos_x, pos_y, pos_z;
    std::vector<double> vel_x, vel_y, vel_z;
    std::vector<int> atom_types;

    // Force field
    LennardJonesParams lj;
    double lj_cutoff{1.2};
    bool use_lj{false};
    std::vector<HarmonicBondParams> bonds;
    std::vector<HarmonicAngleParams> angles;
    std::vector<PeriodicDihedralParams> dihedrals;

    // Electrostatics
    std::string coulomb_type{"cutoff"};
    double coulomb_cutoff{1.2};
    int pme_order{4};
    double pme_grid_spacing{0.12};
    double ewald_coeff{0.34};

    // Output
    std::string trajectory_path{"trajectory.h5"};
    int nstxout{0};
    int nstvout{0};
    int nstenergy{10};
    std::string energy_path{"energy.h5"};
    int checkpoint_interval{500};
    std::string checkpoint_path{"checkpoint.bin"};

    // Coupling
    ThermostatConfig thermostat;
    BarostatConfig barostat;

    // Constraints
    std::string constraint_type{"none"};
    double constraint_tolerance{1e-6};

    // Exclusions (1-2, 1-3, 1-4 pairs to exclude from non-bonded)
    std::vector<int> excl_i;
    std::vector<int> excl_j;
};

class ForceEngine;
class SimulationEngine;

std::unique_ptr<ForceEngine> build_force_engine(const SimulationConfig& cfg);
SimulationEngine build_simulation(const SimulationConfig& cfg);
