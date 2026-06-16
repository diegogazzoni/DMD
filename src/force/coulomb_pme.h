#pragma once

#include "force/force_component.h"
#include <vector>

struct PMEGridParams {
    int nx{32}, ny{32}, nz{32};
    int spline_order{4};
};

class CoulombPME : public ForceComponent {
public:
    explicit CoulombPME(PMEGridParams params);
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
    PMEGridParams params_;
    // FFTW plans and grid buffers will be added when FFTW is integrated
};
