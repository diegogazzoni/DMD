#include "system_data.h"

SystemData::SystemData()
    : n_atoms(0), potential_energy(0.0), step(0), time(0.0) {}

SystemData::SystemData(size_t n)
    : n_atoms(n)
    , pos_x(n, 0.0), pos_y(n, 0.0), pos_z(n, 0.0)
    , vel_x(n, 0.0), vel_y(n, 0.0), vel_z(n, 0.0)
    , forces_x(n, 0.0), forces_y(n, 0.0), forces_z(n, 0.0)
    , masses(n, 0.0), charges(n, 0.0), atom_types(n, 0)
    , potential_energy(0.0), step(0), time(0.0) {}

std::span<const double> SystemData::x() const noexcept { return pos_x; }
std::span<const double> SystemData::y() const noexcept { return pos_y; }
std::span<const double> SystemData::z() const noexcept { return pos_z; }
std::span<double> SystemData::fx() noexcept { return forces_x; }
std::span<double> SystemData::fy() noexcept { return forces_y; }
std::span<double> SystemData::fz() noexcept { return forces_z; }
std::span<double> SystemData::vx() noexcept { return vel_x; }
std::span<double> SystemData::vy() noexcept { return vel_y; }
std::span<double> SystemData::vz() noexcept { return vel_z; }
