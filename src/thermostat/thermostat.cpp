#include "thermostat.h"
#include "core/system_data.h"
#include "core/constants.h"
#include <cmath>

double kinetic_energy(const SystemData& sys) {
    double ekin = 0.0;
    for (size_t i = 0; i < sys.n_atoms; ++i) {
        ekin += 0.5 * sys.masses[i] * (
            sys.vel_x[i] * sys.vel_x[i] +
            sys.vel_y[i] * sys.vel_y[i] +
            sys.vel_z[i] * sys.vel_z[i]);
    }
    return ekin;
}

double temperature(double ekin, size_t n_atoms) {
    if (n_atoms < 2) return 0.0;
    double dof = 3.0 * static_cast<double>(n_atoms) - 3.0;
    return 2.0 * ekin / (dof * dmd_constants::KB);
}

double temperature(const SystemData& sys) {
    return temperature(kinetic_energy(sys), sys.n_atoms);
}
