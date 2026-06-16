#pragma once

#include "force/force_component.h"
#include <span>
#include <vector>

struct HarmonicAngleParams {
    size_t n{0};
    std::vector<int> i, j, k;
    std::vector<double> k_theta, theta0;
};

class HarmonicAngle : public ForceComponent {
public:
    explicit HarmonicAngle(HarmonicAngleParams params);
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
    HarmonicAngleParams params_;
};
