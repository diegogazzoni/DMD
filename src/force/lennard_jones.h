#pragma once

#include "force/force_component.h"
#include <span>
#include <vector>

struct LennardJonesParams {
    size_t n_types{0};
    std::vector<double> sigma;
    std::vector<double> epsilon;
};

class LennardJones : public ForceComponent {
public:
    explicit LennardJones(LennardJonesParams params, double cutoff = 1.2);
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

    void rebuild_pair_list(
        std::span<const double> pos_x,
        std::span<const double> pos_y,
        std::span<const double> pos_z,
        const Cell& cell
    );

private:
    LennardJonesParams params_;
    double cutoff_;
    double cutoff2_;
    std::vector<int> pairs_i_, pairs_j_;
    int rebuild_interval_;
    int step_since_rebuild_;
};
