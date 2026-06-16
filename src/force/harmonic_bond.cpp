#include "harmonic_bond.h"
#include "core/cell.h"
#include <cmath>

HarmonicBond::HarmonicBond(HarmonicBondParams params)
    : params_(std::move(params)) {}

std::string HarmonicBond::name() const {
    return "HarmonicBond";
}

void HarmonicBond::compute(
    std::span<const double> pos_x,
    std::span<const double> pos_y,
    std::span<const double> pos_z,
    const Cell& cell,
    std::span<double> forces_x,
    std::span<double> forces_y,
    std::span<double> forces_z,
    double& energy_out)
{
    double energy = 0.0;
    for (size_t n = 0; n < params_.n; ++n) {
        int ai = params_.i[n];
        int aj = params_.j[n];

        double dx = pos_x[aj] - pos_x[ai];
        double dy = pos_y[aj] - pos_y[ai];
        double dz = pos_z[aj] - pos_z[ai];
        cell.minimum_image(dx, dy, dz, dx, dy, dz);

        double r2 = dx*dx + dy*dy + dz*dz;
        double r = std::sqrt(r2);
        double dr = r - params_.r0[n];
        double f_scalar = -params_.k[n] * dr / r;

        double fx = f_scalar * dx;
        double fy = f_scalar * dy;
        double fz = f_scalar * dz;

        forces_x[ai] -= fx; forces_y[ai] -= fy; forces_z[ai] -= fz;
        forces_x[aj] += fx; forces_y[aj] += fy; forces_z[aj] += fz;

        energy += 0.5 * params_.k[n] * dr * dr;
    }
    energy_out += energy;
}
