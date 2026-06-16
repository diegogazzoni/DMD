#pragma once

#include <span>
#include <string>

struct Cell;

class ForceComponent {
public:
    virtual ~ForceComponent() = default;
    virtual std::string name() const = 0;
    virtual void compute(
        std::span<const double> pos_x,
        std::span<const double> pos_y,
        std::span<const double> pos_z,
        const Cell& cell,
        std::span<double> forces_x,
        std::span<double> forces_y,
        std::span<double> forces_z,
        double& energy_out
    ) = 0;
};
