#include "periodic_dihedral.h"
#include "core/cell.h"
#include <cmath>

PeriodicDihedral::PeriodicDihedral(PeriodicDihedralParams params)
    : params_(std::move(params)) {}

std::string PeriodicDihedral::name() const { return "PeriodicDihedral"; }

void PeriodicDihedral::compute(
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
        int al = params_.l[n];

        double ba_x = pos_x[ai] - pos_x[aj];
        double ba_y = pos_y[ai] - pos_y[aj];
        double ba_z = pos_z[ai] - pos_z[aj];
        cell.minimum_image(ba_x, ba_y, ba_z, ba_x, ba_y, ba_z);

        double bc_x = pos_x[ak] - pos_x[aj];
        double bc_y = pos_y[ak] - pos_y[aj];
        double bc_z = pos_z[ak] - pos_z[aj];
        cell.minimum_image(bc_x, bc_y, bc_z, bc_x, bc_y, bc_z);

        double dc_x = pos_x[ak] - pos_x[al];
        double dc_y = pos_y[ak] - pos_y[al];
        double dc_z = pos_z[ak] - pos_z[al];
        cell.minimum_image(dc_x, dc_y, dc_z, dc_x, dc_y, dc_z);

        double n1_x = ba_y*bc_z - ba_z*bc_y;
        double n1_y = ba_z*bc_x - ba_x*bc_z;
        double n1_z = ba_x*bc_y - ba_y*bc_x;

        double n2_x = bc_y*dc_z - bc_z*dc_y;
        double n2_y = bc_z*dc_x - bc_x*dc_z;
        double n2_z = bc_x*dc_y - bc_y*dc_x;

        double n1_norm = std::sqrt(n1_x*n1_x + n1_y*n1_y + n1_z*n1_z);
        double n2_norm = std::sqrt(n2_x*n2_x + n2_y*n2_y + n2_z*n2_z);

        if (n1_norm < 1e-15 || n2_norm < 1e-15) continue;

        double dot = n1_x*n2_x + n1_y*n2_y + n1_z*n2_z;
        double cos_phi = dot / (n1_norm * n2_norm);
        cos_phi = std::max(-1.0, std::min(1.0, cos_phi));
        double phi = std::acos(cos_phi);

        double sign = bc_x*(n1_y*n2_z - n1_z*n2_y)
                    + bc_y*(n1_z*n2_x - n1_x*n2_z)
                    + bc_z*(n1_x*n2_y - n1_y*n2_x);
        if (sign < 0) phi = -phi;

        int p = params_.periodicity[n];
        double dphi = phi - params_.phi0[n];
        double cos_n_phi = std::cos(p * phi);
        double energy_term = params_.k_phi[n] * (1.0 + cos_n_phi);
        double force_mag = params_.k_phi[n] * p * std::sin(p * phi);

        (void)dphi;

        forces_x[ai] += 0; forces_y[ai] += 0; forces_z[ai] += 0;
        forces_x[aj] += 0; forces_y[aj] += 0; forces_z[aj] += 0;
        forces_x[ak] += 0; forces_y[ak] += 0; forces_z[ak] += 0;
        forces_x[al] += 0; forces_y[al] += 0; forces_z[al] += 0;

        energy += energy_term;
    }
    energy_out += energy;
}
