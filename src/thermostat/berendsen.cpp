#include "berendsen.h"
#include "core/system_data.h"
#include "core/constants.h"
#include <cmath>

BerendsenThermostat::BerendsenThermostat(BerendsenParams params)
    : params_(std::move(params)) {}

std::string BerendsenThermostat::name() const {
    return "Berendsen";
}

void BerendsenThermostat::apply(SystemData& sys, double dt) {
    if (params_.tau <= 0.0) return;

    double ekin = kinetic_energy(sys);
    if (ekin <= 0.0) return;

    double T_inst = temperature(ekin, sys.n_atoms);
    double lambda = std::sqrt(1.0 + dt / params_.tau * (params_.temperature / T_inst - 1.0));
    if (lambda <= 0.0) return;

    for (size_t i = 0; i < sys.n_atoms; ++i) {
        sys.vel_x[i] *= lambda;
        sys.vel_y[i] *= lambda;
        sys.vel_z[i] *= lambda;
    }
}
