#pragma once

#include "force/force_component.h"
#include <vector>

class CoulombExclusion : public ForceComponent {
public:
    std::string name() const override;
    void compute(
        std::span<const double> pos_x,
        std::span<const double> pos_y,
        std::span<const double> pos_z,
        const Cell& cell,
        std::span<double> forces_x,
        std::span<double> forces_y,
        std::span<double> forces_z,
        double& energy_out
    ) override;

    void set_charges(std::span<const double> charges);
    void set_exclusions(std::span<const int> excl_i, std::span<const int> excl_j);

private:
    std::vector<double> charges_;
    std::vector<int> excl_i_, excl_j_;
};
