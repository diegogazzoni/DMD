#include "barostat.h"
#include "core/system_data.h"
#include "core/cell.h"
#include <cmath>

double virial(const SystemData& sys) {
    double w = 0.0;
    for (size_t i = 0; i < sys.n_atoms; ++i) {
        w += sys.pos_x[i] * sys.forces_x[i] +
             sys.pos_y[i] * sys.forces_y[i] +
             sys.pos_z[i] * sys.forces_z[i];
    }
    return w;
}

double pressure(double ekin, double virial_val, double volume) {
    if (volume <= 0.0) return 0.0;
    return (2.0 * ekin + virial_val) / (3.0 * volume);
}

double pressure(const SystemData& sys, double volume) {
    double ekin = 0.0;
    for (size_t i = 0; i < sys.n_atoms; ++i) {
        ekin += 0.5 * sys.masses[i] * (
            sys.vel_x[i] * sys.vel_x[i] +
            sys.vel_y[i] * sys.vel_y[i] +
            sys.vel_z[i] * sys.vel_z[i]);
    }
    return pressure(ekin, virial(sys), volume);
}
