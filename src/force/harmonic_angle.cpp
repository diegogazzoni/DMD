#include "harmonic_angle.h"
#include "core/cell.h"
#include <cmath>

HarmonicAngle::HarmonicAngle(HarmonicAngleParams params)
    : params_(std::move(params)) {}

std::string HarmonicAngle::name() const { return "HarmonicAngle"; }

void HarmonicAngle::compute(
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
        int ak = params_.k[n];

        double dx_ij = pos_x[ai] - pos_x[aj];
        double dy_ij = pos_y[ai] - pos_y[aj];
        double dz_ij = pos_z[ai] - pos_z[aj];
        cell.minimum_image(dx_ij, dy_ij, dz_ij, dx_ij, dy_ij, dz_ij);

        double dx_kj = pos_x[ak] - pos_x[aj];
        double dy_kj = pos_y[ak] - pos_y[aj];
        double dz_kj = pos_z[ak] - pos_z[aj];
        cell.minimum_image(dx_kj, dy_kj, dz_kj, dx_kj, dy_kj, dz_kj);

        double r_ij2 = dx_ij*dx_ij + dy_ij*dy_ij + dz_ij*dz_ij;
        double r_kj2 = dx_kj*dx_kj + dy_kj*dy_kj + dz_kj*dz_kj;
        double r_ij = std::sqrt(r_ij2);
        double r_kj = std::sqrt(r_kj2);

        double dot = dx_ij*dx_kj + dy_ij*dy_kj + dz_ij*dz_kj;
        double cos_theta = dot / (r_ij * r_kj);
        cos_theta = std::max(-1.0, std::min(1.0, cos_theta));
        double theta = std::acos(cos_theta);
        double dtheta = theta - params_.theta0[n];
        double f_scalar = -params_.k_theta[n] * dtheta / std::sqrt(1.0 - cos_theta*cos_theta);

        double fx_ij = f_scalar * (dx_kj / (r_ij * r_kj) - cos_theta * dx_ij / r_ij2);
        double fy_ij = f_scalar * (dy_kj / (r_ij * r_kj) - cos_theta * dy_ij / r_ij2);
        double fz_ij = f_scalar * (dz_kj / (r_ij * r_kj) - cos_theta * dz_ij / r_ij2);

        double fx_kj = f_scalar * (dx_ij / (r_ij * r_kj) - cos_theta * dx_kj / r_kj2);
        double fy_kj = f_scalar * (dy_ij / (r_ij * r_kj) - cos_theta * dy_kj / r_kj2);
        double fz_kj = f_scalar * (dz_ij / (r_ij * r_kj) - cos_theta * dz_kj / r_kj2);

        forces_x[ai] += fx_ij; forces_y[ai] += fy_ij; forces_z[ai] += fz_ij;
        forces_x[aj] -= fx_ij + fx_kj; forces_y[aj] -= fy_ij + fy_kj; forces_z[aj] -= fz_ij + fz_kj;
        forces_x[ak] += fx_kj; forces_y[ak] += fy_kj; forces_z[ak] += fz_kj;

        energy += 0.5 * params_.k_theta[n] * dtheta * dtheta;
    }
    energy_out += energy;
}
