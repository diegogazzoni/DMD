#pragma once

#include <cstddef>
#include <span>
#include <vector>

struct SystemData {
    size_t n_atoms;

    std::vector<double> pos_x, pos_y, pos_z;
    std::vector<double> vel_x, vel_y, vel_z;
    std::vector<double> forces_x, forces_y, forces_z;
    std::vector<double> masses;
    std::vector<double> charges;
    std::vector<int>    atom_types;

    double potential_energy;
    size_t step;
    double time;

    SystemData();
    explicit SystemData(size_t n);

    std::span<const double> x() const noexcept;
    std::span<const double> y() const noexcept;
    std::span<const double> z() const noexcept;
    std::span<double> fx() noexcept;
    std::span<double> fy() noexcept;
    std::span<double> fz() noexcept;
    std::span<double> vx() noexcept;
    std::span<double> vy() noexcept;
    std::span<double> vz() noexcept;
};
