#include "integrator.h"
#include "core/system_data.h"

void Integrator::step(SystemData& sys, double dt) {
    size_t n = sys.n_atoms;
    for (size_t i = 0; i < n; ++i) {
        double inv_mass = 1.0 / sys.masses[i];
        sys.vel_x[i] += 0.5 * sys.forces_x[i] * inv_mass * dt;
        sys.vel_y[i] += 0.5 * sys.forces_y[i] * inv_mass * dt;
        sys.vel_z[i] += 0.5 * sys.forces_z[i] * inv_mass * dt;

        sys.pos_x[i] += sys.vel_x[i] * dt;
        sys.pos_y[i] += sys.vel_y[i] * dt;
        sys.pos_z[i] += sys.vel_z[i] * dt;
    }

    for (size_t i = 0; i < n; ++i) {
        double inv_mass = 1.0 / sys.masses[i];
        sys.vel_x[i] += 0.5 * sys.forces_x[i] * inv_mass * dt;
        sys.vel_y[i] += 0.5 * sys.forces_y[i] * inv_mass * dt;
        sys.vel_z[i] += 0.5 * sys.forces_z[i] * inv_mass * dt;
    }
}
