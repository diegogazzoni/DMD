#pragma once

#include "core/cell.h"
#include "force/lennard_jones.h"
#include "force/harmonic_bond.h"
#include "force/harmonic_angle.h"
#include "force/periodic_dihedral.h"
#include <memory>
#include <string>
#include <vector>

struct SimulationConfig {
    double dt{0.002};
    int n_steps{1000};

    double box_size{3.0};

    std::vector<double> masses;
    std::vector<double> charges;
    std::vector<double> pos_x, pos_y, pos_z;
    std::vector<double> vel_x, vel_y, vel_z;
    std::vector<int> atom_types;

    LennardJonesParams lj;
    double lj_cutoff{1.2};
    bool use_lj{false};

    std::vector<HarmonicBondParams> bonds;
    std::vector<HarmonicAngleParams> angles;
    std::vector<PeriodicDihedralParams> dihedrals;

    int trajectory_interval{100};
    int checkpoint_interval{500};
    std::string checkpoint_path{"checkpoint.bin"};
};

class ForceEngine;
class SimulationEngine;

std::unique_ptr<ForceEngine> build_force_engine(const SimulationConfig& cfg);
SimulationEngine build_simulation(const SimulationConfig& cfg);
