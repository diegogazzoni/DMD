#include "lennard_jones.h"
#include "core/cell.h"
#include <cmath>

LennardJones::LennardJones(LennardJonesParams params, double cutoff)
    : params_(std::move(params))
    , cutoff_(cutoff)
    , cutoff2_(cutoff * cutoff)
    , rebuild_interval_(10)
    , step_since_rebuild_(0) {}

void LennardJones::set_atom_types(std::span<const int> types) {
    atom_types_.assign(types.begin(), types.end());
}

std::string LennardJones::name() const { return "LennardJones"; }

void LennardJones::rebuild_pair_list(
    std::span<const double> pos_x,
    std::span<const double> pos_y,
    std::span<const double> pos_z,
    const Cell& cell)
{
    pairs_i_.clear();
    pairs_j_.clear();
    size_t n = pos_x.size();
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            double dx = pos_x[j] - pos_x[i];
            double dy = pos_y[j] - pos_y[i];
            double dz = pos_z[j] - pos_z[i];
            cell.minimum_image(dx, dy, dz, dx, dy, dz);
            double r2 = dx*dx + dy*dy + dz*dz;
            if (r2 < cutoff2_) {
                pairs_i_.push_back(static_cast<int>(i));
                pairs_j_.push_back(static_cast<int>(j));
            }
        }
    }
    step_since_rebuild_ = 0;
}

void LennardJones::compute(
    std::span<const double> pos_x,
    std::span<const double> pos_y,
    std::span<const double> pos_z,
    const Cell& cell,
    std::span<double> forces_x,
    std::span<double> forces_y,
    std::span<double> forces_z,
    double& energy_out)
{
    if (step_since_rebuild_ >= rebuild_interval_) {
        rebuild_pair_list(pos_x, pos_y, pos_z, cell);
    }

    double energy = 0.0;
    for (size_t p = 0; p < pairs_i_.size(); ++p) {
        int i = pairs_i_[p];
        int j = pairs_j_[p];

        double dx = pos_x[j] - pos_x[i];
        double dy = pos_y[j] - pos_y[i];
        double dz = pos_z[j] - pos_z[i];
        cell.minimum_image(dx, dy, dz, dx, dy, dz);

        double r2 = dx*dx + dy*dy + dz*dz;
        if (r2 >= cutoff2_) continue;

        int ti = atom_types_[i];
        int tj = atom_types_[j];
        double eps = params_.epsilon[ti * params_.n_types + tj];
        double sig = params_.sigma[ti * params_.n_types + tj];

        double sig2 = sig * sig;
        double sig6 = sig2 * sig2 * sig2;
        double r6 = r2 * r2 * r2;
        double sr6 = sig6 / r6;
        double f_scalar = 24.0 * eps * (2.0 * sr6 * sr6 - sr6) / std::sqrt(r2);

        forces_x[i] += f_scalar * dx;
        forces_y[i] += f_scalar * dy;
        forces_z[i] += f_scalar * dz;
        forces_x[j] -= f_scalar * dx;
        forces_y[j] -= f_scalar * dy;
        forces_z[j] -= f_scalar * dz;

        energy += 4.0 * eps * (sr6 * sr6 - sr6);
    }
    energy_out += energy;
    ++step_since_rebuild_;
}
