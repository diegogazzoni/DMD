#include "coulomb_direct.h"
#include "core/cell.h"
#include <cmath>

CoulombDirect::CoulombDirect(CoulombDirectParams params)
    : params_(params), cutoff2_(params.cutoff * params.cutoff) {}

std::string CoulombDirect::name() const { return "CoulombDirect"; }

void CoulombDirect::set_charges(std::span<const double> charges) {
    charges_.assign(charges.begin(), charges.end());
}

void CoulombDirect::rebuild_pair_list(
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
}

void CoulombDirect::compute(
    std::span<const double> pos_x,
    std::span<const double> pos_y,
    std::span<const double> pos_z,
    const Cell& cell,
    std::span<double> forces_x,
    std::span<double> forces_y,
    std::span<double> forces_z,
    double& energy_out)
{
#if defined(__APPLE__)
    double k_e = 138.935456; // kJ/mol * nm / e^2
#else
    double k_e = 138.935456;
#endif
    double alpha = params_.ewald_coefficient;
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

        double r = std::sqrt(r2);
        double qi = charges_[i];
        double qj = charges_[j];
        double qiqj = qi * qj;

        double erfc_ar = std::erfc(alpha * r);
        double force_scalar = k_e * qiqj * (erfc_ar / (r2 * r) + 2.0 * alpha * std::exp(-alpha * alpha * r2) / (std::sqrt(M_PI) * r2));

        forces_x[i] += force_scalar * dx;
        forces_y[i] += force_scalar * dy;
        forces_z[i] += force_scalar * dz;
        forces_x[j] -= force_scalar * dx;
        forces_y[j] -= force_scalar * dy;
        forces_z[j] -= force_scalar * dz;

        energy += k_e * qiqj * erfc_ar / r;
    }
    energy_out += energy;
}
