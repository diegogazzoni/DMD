#pragma once

#include "force/force_component.h"
#include <vector>

struct CoulombDirectParams {
    double cutoff{1.2};
    double ewald_coefficient{3.5};
};

class CoulombDirect : public ForceComponent {
public:
    explicit CoulombDirect(CoulombDirectParams params);
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
    void rebuild_pair_list(
        std::span<const double> pos_x,
        std::span<const double> pos_y,
        std::span<const double> pos_z,
        const Cell& cell
    );

private:
    CoulombDirectParams params_;
    std::vector<int> pairs_i_, pairs_j_;
    std::vector<double> charges_;
    double cutoff2_;
};
