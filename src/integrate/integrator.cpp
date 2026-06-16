#include "integrator.h"
#include "core/system_data.h"

void Integrator::half_kick(SystemData& sys, double dt) {
    double half_dt = 0.5 * dt;
    for (size_t i = 0; i < sys.n_atoms; ++i) {
        double inv_mass = 1.0 / sys.masses[i];
        sys.vel_x[i] += sys.forces_x[i] * inv_mass * half_dt;
        sys.vel_y[i] += sys.forces_y[i] * inv_mass * half_dt;
        sys.vel_z[i] += sys.forces_z[i] * inv_mass * half_dt;
    }
}

void Integrator::advance(SystemData& sys, double dt) {
    for (size_t i = 0; i < sys.n_atoms; ++i) {
        sys.pos_x[i] += sys.vel_x[i] * dt;
        sys.pos_y[i] += sys.vel_y[i] * dt;
        sys.pos_z[i] += sys.vel_z[i] * dt;
    }
}
