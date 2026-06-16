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

        // ba = r_i - r_j, bc = r_k - r_j, dc = r_k - r_l
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

        // N1 = ba × bc, N2 = bc × dc
        double n1_x = ba_y*bc_z - ba_z*bc_y;
        double n1_y = ba_z*bc_x - ba_x*bc_z;
        double n1_z = ba_x*bc_y - ba_y*bc_x;

        double n2_x = bc_y*dc_z - bc_z*dc_y;
        double n2_y = bc_z*dc_x - bc_x*dc_z;
        double n2_z = bc_x*dc_y - bc_y*dc_x;

        double n1_sq = n1_x*n1_x + n1_y*n1_y + n1_z*n1_z;
        double n2_sq = n2_x*n2_x + n2_y*n2_y + n2_z*n2_z;
        double bc_sq = bc_x*bc_x + bc_y*bc_y + bc_z*bc_z;

        if (n1_sq < 1e-30 || n2_sq < 1e-30 || bc_sq < 1e-30) continue;

        double n1_norm = std::sqrt(n1_sq);
        double n2_norm = std::sqrt(n2_sq);
        double bc_norm = std::sqrt(bc_sq);

        // cos_phi and sin_phi
        double dot = n1_x*n2_x + n1_y*n2_y + n1_z*n2_z;
        double cos_phi = dot / (n1_norm * n2_norm);
        cos_phi = std::max(-1.0, std::min(1.0, cos_phi));

        double triple = bc_x*(n1_y*n2_z - n1_z*n2_y)
                      + bc_y*(n1_z*n2_x - n1_x*n2_z)
                      + bc_z*(n1_x*n2_y - n1_y*n2_x);
        double sin_phi = triple / (bc_norm * n1_norm * n2_norm);
        sin_phi = std::max(-1.0, std::min(1.0, sin_phi));

        double phi = std::atan2(sin_phi, cos_phi);

        int p = params_.periodicity[n];
        double k = params_.k_phi[n];
        double phi0 = params_.phi0[n];

        double cos_np = std::cos(p * phi - phi0);
        energy += k * (1.0 + cos_np);

        // Force prefactor: -dE/dphi = k * p * sin(p*phi - phi0)
        double pref = k * p * std::sin(p * phi - phi0);
        if (std::abs(pref) < 1e-30) continue;

        // Intermediate vectors for gradient computation
        // dphi/dr_i = (bc × N1) / n1_sq
        double gi_x = (bc_y*n1_z - bc_z*n1_y) / n1_sq;
        double gi_y = (bc_z*n1_x - bc_x*n1_z) / n1_sq;
        double gi_z = (bc_x*n1_y - bc_y*n1_x) / n1_sq;

        // dphi/dr_l = -(bc × N2) / n2_sq
        double gl_x = -(bc_y*n2_z - bc_z*n2_y) / n2_sq;
        double gl_y = -(bc_z*n2_x - bc_x*n2_z) / n2_sq;
        double gl_z = -(bc_x*n2_y - bc_y*n2_x) / n2_sq;

        // dphi/dr_j = (N1 × ba)/n1_sq + (dc × N2)/n2_sq - (bc × N1)/n1_sq
        double n1xba_x = n1_y*ba_z - n1_z*ba_y;
        double n1xba_y = n1_z*ba_x - n1_x*ba_z;
        double n1xba_z = n1_x*ba_y - n1_y*ba_x;

        double dcxn2_x = dc_y*n2_z - dc_z*n2_y;
        double dcxn2_y = dc_z*n2_x - dc_x*n2_z;
        double dcxn2_z = dc_x*n2_y - dc_y*n2_x;

        double gj_x = n1xba_x / n1_sq + dcxn2_x / n2_sq - gi_x;
        double gj_y = n1xba_y / n1_sq + dcxn2_y / n2_sq - gi_y;
        double gj_z = n1xba_z / n1_sq + dcxn2_z / n2_sq - gi_z;

        // dphi/dr_k = (bc × N2)/n2_sq - (N1 × ba)/n1_sq - (dc × N2)/n2_sq
        double bxn2_x = (bc_y*n2_z - bc_z*n2_y) / n2_sq;
        double bxn2_y = (bc_z*n2_x - bc_x*n2_z) / n2_sq;
        double bxn2_z = (bc_x*n2_y - bc_y*n2_x) / n2_sq;

        double gk_x = bxn2_x - n1xba_x / n1_sq - dcxn2_x / n2_sq;
        double gk_y = bxn2_y - n1xba_y / n1_sq - dcxn2_y / n2_sq;
        double gk_z = bxn2_z - n1xba_z / n1_sq - dcxn2_z / n2_sq;

        // Forces: F_m = pref * dphi/dr_m
        forces_x[ai] += pref * gi_x;
        forces_y[ai] += pref * gi_y;
        forces_z[ai] += pref * gi_z;

        forces_x[aj] += pref * gj_x;
        forces_y[aj] += pref * gj_y;
        forces_z[aj] += pref * gj_z;

        forces_x[ak] += pref * gk_x;
        forces_y[ak] += pref * gk_y;
        forces_z[ak] += pref * gk_z;

        forces_x[al] += pref * gl_x;
        forces_y[al] += pref * gl_y;
        forces_z[al] += pref * gl_z;
    }
    energy_out += energy;
}
