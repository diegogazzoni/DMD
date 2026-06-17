#include "coulomb_exclusion.h"
#include "core/cell.h"
#include <cmath>

std::string CoulombExclusion::name() const { return "CoulombExclusion"; }

void CoulombExclusion::set_charges(std::span<const double> charges) {
    charges_.assign(charges.begin(), charges.end());
}

void CoulombExclusion::set_exclusions(std::span<const int> excl_i, std::span<const int> excl_j) {
    excl_i_.assign(excl_i.begin(), excl_i.end());
    excl_j_.assign(excl_j.begin(), excl_j.end());
}

void CoulombExclusion::compute(
    std::span<const double> pos_x,
    std::span<const double> pos_y,
    std::span<const double> pos_z,
    const Cell& cell,
    std::span<double> forces_x,
    std::span<double> forces_y,
    std::span<double> forces_z,
    double& energy_out)
{
    double k_e = 138.935456;
    double energy = 0.0;

    for (size_t p = 0; p < excl_i_.size(); ++p) {
        int i = excl_i_[p];
        int j = excl_j_[p];

        double dx = pos_x[j] - pos_x[i];
        double dy = pos_y[j] - pos_y[i];
        double dz = pos_z[j] - pos_z[i];
        cell.minimum_image(dx, dy, dz, dx, dy, dz);

        double r2 = dx*dx + dy*dy + dz*dz;
        double r = std::sqrt(r2);

        double qi = charges_[i];
        double qj = charges_[j];
        double qiqj = qi * qj;

        double force_scalar = k_e * qiqj / (r2 * r);
        forces_x[i] -= force_scalar * dx;
        forces_y[i] -= force_scalar * dy;
        forces_z[i] -= force_scalar * dz;
        forces_x[j] += force_scalar * dx;
        forces_y[j] += force_scalar * dy;
        forces_z[j] += force_scalar * dz;

        energy += k_e * qiqj / r;
    }
    energy_out -= energy;
}
