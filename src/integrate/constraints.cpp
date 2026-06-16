#include "constraints.h"
#include "core/system_data.h"
#include <cmath>

Constraints::Constraints(ConstraintsParams params)
    : params_(std::move(params)) {}

void Constraints::apply(SystemData& sys, double dt) {
    // SHAKE: iterative constraint correction (post-position)
    const double tolerance = 1e-8;
    const int max_iterations = 100;

    for (int iter = 0; iter < max_iterations; ++iter) {
        double max_error = 0.0;
        for (size_t n = 0; n < params_.n; ++n) {
            int ai = params_.i[n];
            int aj = params_.j[n];

            double dx = sys.pos_x[aj] - sys.pos_x[ai];
            double dy = sys.pos_y[aj] - sys.pos_y[ai];
            double dz = sys.pos_z[aj] - sys.pos_z[ai];

            double r2 = dx*dx + dy*dy + dz*dz;
            double d2 = params_.distance_target[n] * params_.distance_target[n];
            double error = r2 - d2;

            if (std::abs(error) > max_error) max_error = std::abs(error);

            if (std::abs(error) > tolerance) {
                double inv_mass_i = 1.0 / sys.masses[ai];
                double inv_mass_j = 1.0 / sys.masses[aj];
                double inv_mass_sum = inv_mass_i + inv_mass_j;
                double lambda = error / (2.0 * d2 * inv_mass_sum);

                double fx = lambda * dx;
                double fy = lambda * dy;
                double fz = lambda * dz;

                sys.pos_x[ai] += fx * inv_mass_i;
                sys.pos_y[ai] += fy * inv_mass_i;
                sys.pos_z[ai] += fz * inv_mass_i;
                sys.pos_x[aj] -= fx * inv_mass_j;
                sys.pos_y[aj] -= fy * inv_mass_j;
                sys.pos_z[aj] -= fz * inv_mass_j;
            }
        }
        if (max_error < tolerance) break;
    }
}
