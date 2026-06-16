#include "andersen_barostat.h"
#include "core/system_data.h"
#include "core/cell.h"
#include "core/constants.h"
#include "thermostat/thermostat.h"
#include <cmath>

AndersenBarostat::AndersenBarostat(AndersenBarostatParams params, size_t n_atoms)
    : params_(std::move(params))
{
    double dof = 3.0 * static_cast<double>(n_atoms) - 3.0;
    if (dof < 1.0) dof = 1.0;
    W_ = dof * dmd_constants::KB * params_.pressure * params_.tau * params_.tau;
}

std::string AndersenBarostat::name() const {
    return "AndersenBarostat";
}

void AndersenBarostat::apply(SystemData& sys, Cell& cell, double dt) {
    double volume = cell.volume();
    if (volume <= 0.0) return;

    double ekin = kinetic_energy(sys);
    double P_inst = pressure(ekin, virial(sys), volume);

    // update piston velocity: eps_dot += dt/2 * 3V*(P-P0)/W
    double forcing = 3.0 * volume * (P_inst - params_.pressure) / W_;
    eps_dot_ = eps_dot_ + 0.5 * dt * forcing;

    // scale positions by exp(eps_dot * dt/2)
    // This is the "advance" step for the barostat
    // In the full MTK scheme, this is split across the MD step
    double scale = std::exp(eps_dot_ * dt);
    for (size_t i = 0; i < sys.n_atoms; ++i) {
        sys.pos_x[i] *= scale;
        sys.pos_y[i] *= scale;
        sys.pos_z[i] *= scale;
    }

    cell.scale(scale, scale, scale);

    // complete velocity update (second half)
    eps_dot_ = eps_dot_ + 0.5 * dt * forcing;
}
