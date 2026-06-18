#include "constraints.h"
#include "core/cell.h"
#include "core/system_data.h"
#include <cmath>

Constraints::Constraints(ConstraintsParams params)
    : params_(std::move(params)) {}

void Constraints::apply(SystemData& sys, const Cell& cell) {
    const int max_iterations = 100;
    const double tolerance = params_.tolerance;
    size_t n = params_.i.size();

    for (int iter = 0; iter < max_iterations; ++iter) {
        double max_error = 0.0;
        for (size_t c = 0; c < n; ++c) {
            int ai = params_.i[c];
            int aj = params_.j[c];

            double dx = sys.pos_x[aj] - sys.pos_x[ai];
            double dy = sys.pos_y[aj] - sys.pos_y[ai];
            double dz = sys.pos_z[aj] - sys.pos_z[ai];
            cell.minimum_image(dx, dy, dz, dx, dy, dz);

            double r2 = dx*dx + dy*dy + dz*dz;
            double d2 = params_.distance_target[c] * params_.distance_target[c];
            double error = r2 - d2;

            if (std::abs(error) > max_error) max_error = std::abs(error);

            if (std::abs(error) > tolerance) {
                double inv_mass_i = 1.0 / sys.masses[ai];
                double inv_mass_j = 1.0 / sys.masses[aj];
                double inv_mass_sum = inv_mass_i + inv_mass_j;
                double lambda = error / (2.0 * d2 * inv_mass_sum);

                sys.pos_x[ai] += lambda * dx * inv_mass_i;
                sys.pos_y[ai] += lambda * dy * inv_mass_i;
                sys.pos_z[ai] += lambda * dz * inv_mass_i;
                sys.pos_x[aj] -= lambda * dx * inv_mass_j;
                sys.pos_y[aj] -= lambda * dy * inv_mass_j;
                sys.pos_z[aj] -= lambda * dz * inv_mass_j;
            }
        }
        if (max_error < tolerance) break;
    }
}

void Constraints::apply_rattle(SystemData& sys, const Cell& cell) {
    const int max_iterations = 100;
    const double tolerance = params_.tolerance;
    size_t n = params_.i.size();

    for (int iter = 0; iter < max_iterations; ++iter) {
        double max_error = 0.0;
        for (size_t c = 0; c < n; ++c) {
            int ai = params_.i[c];
            int aj = params_.j[c];

            double dx = sys.pos_x[aj] - sys.pos_x[ai];
            double dy = sys.pos_y[aj] - sys.pos_y[ai];
            double dz = sys.pos_z[aj] - sys.pos_z[ai];
            cell.minimum_image(dx, dy, dz, dx, dy, dz);

            double dvx = sys.vel_x[aj] - sys.vel_x[ai];
            double dvy = sys.vel_y[aj] - sys.vel_y[ai];
            double dvz = sys.vel_z[aj] - sys.vel_z[ai];

            double r_dot_v = dx*dvx + dy*dvy + dz*dvz;
            double r2 = dx*dx + dy*dy + dz*dz;

            double err = std::abs(r_dot_v);
            if (err > max_error) max_error = err;

            if (err > tolerance) {
                double inv_mass_i = 1.0 / sys.masses[ai];
                double inv_mass_j = 1.0 / sys.masses[aj];
                double inv_mass_sum = inv_mass_i + inv_mass_j;
                double lambda = r_dot_v / (r2 * inv_mass_sum);

                sys.vel_x[ai] += lambda * dx * inv_mass_i;
                sys.vel_y[ai] += lambda * dy * inv_mass_i;
                sys.vel_z[ai] += lambda * dz * inv_mass_i;
                sys.vel_x[aj] -= lambda * dx * inv_mass_j;
                sys.vel_y[aj] -= lambda * dy * inv_mass_j;
                sys.vel_z[aj] -= lambda * dz * inv_mass_j;
            }
        }
        if (max_error < tolerance) break;
    }
}
