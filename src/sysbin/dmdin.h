#pragma once

#include "sim/simulation_config.h"
#include <string>

// dmdin binary format for system input
// Magic: "DMDN" (0x444D444E), Version 1
//
// Layout:
//   HEADER (120 bytes) — see DMDinHeader
//   masses[n_types] (double)
//   charges[n_types] (double)
//   FF section at offset_ff (see below)
//   Positions section at offset_pos (see below)
//
// FF section:
//   n_bonds (uint32_t)
//   bonds_i[n_bonds] (int32_t), bonds_j[n_bonds] (int32_t)
//   bonds_k[n_bonds] (double), bonds_r0[n_bonds] (double)
//   n_angles (uint32_t)
//   angles_i[n_angles] (int32_t), angles_j[n_angles] (int32_t), angles_k[n_angles] (int32_t)
//   angles_k_theta[n_angles] (double), angles_theta0[n_angles] (double)
//   n_dihedrals (uint32_t)
//   dih_i[n_dihedrals] (int32_t), dih_j[n_dihedrals] (int32_t),
//   dih_k[n_dihedrals] (int32_t), dih_l[n_dihedrals] (int32_t)
//   dih_per[n_dihedrals] (int32_t)
//   dih_k_phi[n_dihedrals] (double), dih_phi0[n_dihedrals] (double)
//   n_lj_types (uint32_t)
//   sigma[n_lj_types * n_lj_types] (double)
//   epsilon[n_lj_types * n_lj_types] (double)
//
// Positions section:
//   pos_x[n_atoms] (float), pos_y[n_atoms] (float), pos_z[n_atoms] (float)
//   vel_x[n_atoms] (float), vel_y[n_atoms] (float), vel_z[n_atoms] (float) [if pos_flags & 1]
//   atom_types[n_atoms] (int32_t)

SimulationConfig read_dmdin(const std::string& path);
void write_dmdin(const std::string& path, const SimulationConfig& cfg);
