#pragma once

#include "force/force_component.h"
#include <span>
#include <vector>

struct ImproperParams {
    size_t n{0};
    std::vector<int> i, j, k, l;
    std::vector<double> k_phi, phi0;
};

class HarmonicImproper : public ForceComponent {
public:
    explicit HarmonicImproper(ImproperParams params);
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

private:
    ImproperParams params_;
};
