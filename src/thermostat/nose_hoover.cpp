#include "nose_hoover.h"
#include "core/system_data.h"
#include "core/constants.h"
#include <cmath>

NoseHooverThermostat::NoseHooverThermostat(NoseHooverParams params, size_t n_atoms)
    : params_(std::move(params))
{
    double dof = 3.0 * static_cast<double>(n_atoms) - 3.0;
    if (dof < 1.0) dof = 1.0;
    Q_ = dof * dmd_constants::KB * params_.temperature * params_.tau * params_.tau;
}

std::string NoseHooverThermostat::name() const {
    return "NoseHoover";
}

void NoseHooverThermostat::apply(SystemData& sys, double dt) {
    double ekin = kinetic_energy(sys);
    double T_inst = temperature(ekin, sys.n_atoms);

    // update friction: zeta += dt/2 * (T - T0) / Q  (first half)
    zeta_ += 0.5 * dt * dmd_constants::KB * (T_inst - params_.temperature) / Q_;

    // scale velocities: v *= exp(-zeta * dt/2)
    double scale = std::exp(-zeta_ * 0.5 * dt);
    for (size_t i = 0; i < sys.n_atoms; ++i) {
        sys.vel_x[i] *= scale;
        sys.vel_y[i] *= scale;
        sys.vel_z[i] *= scale;
    }
}
