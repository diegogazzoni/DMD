#include "berendsen_barostat.h"
#include "core/system_data.h"
#include "core/cell.h"
#include "thermostat/thermostat.h"
#include <cmath>

BerendsenBarostat::BerendsenBarostat(BerendsenBarostatParams params)
    : params_(std::move(params)) {}

std::string BerendsenBarostat::name() const {
    return "BerendsenBarostat";
}

void BerendsenBarostat::apply(SystemData& sys, Cell& cell, double dt) {
    if (params_.tau <= 0.0) return;

    double volume = cell.volume();
    if (volume <= 0.0) return;

    double ekin = kinetic_energy(sys);
    double P_inst = pressure(ekin, virial(sys), volume);

    // convert bar to kJ/(nm^3) = 100 kJ/(mol*nm^3) approximately
    // 1 bar = 0.0001 kJ/(mol*nm^3) with consistent units
    // In MD units: P(bar) * 0.0001 = P(kJ/(mol*nm^3))
    // But the compressibility also needs consistent units
    // For simplicity, work in internal units and convert:
    // 1 bar = 6.022e-5 kJ/(mol*A^3) approximately
    // Let's use a simpler approach: scale factor directly

    // Berendsen: mu = (1 + dt/tau_P * beta * (P - P0))^(1/3)
    // beta = compressibility in consistent units
    double beta = params_.compressibility;

    // P_inst and params_.pressure are both in bar units from pressure()
    // The conversion factor: 1 bar in internal units (kJ/(mol*nm^3))
    // P(internal) = P(bar) * 0.0602214086 (approx for kJ/(mol*nm^3))
    // The full expression is mu^(1/3) = 1 - gamma*dt/tau_P * (P0 - P)
    // gamma = beta * compressibility factor

    // Simplest: use the standard Berendsen barostat formula
    double mu = 1.0 - (dt / params_.tau) * beta * (params_.pressure - P_inst);

    // For water at 300K, beta ~ 4.5e-5 bar^-1
    // Need to convert to consistent units: 
    // compressibility in bar^-1 * P(bar) is dimensionless
    // The factor converts bar to kJ/(mol*nm^3):
    // P(kJ/(mol*nm^3)) = P(bar) * 0.0602214086
    // So beta(internal) = beta(bar^-1) / 0.0602214086

    // Actually, let's just be pragmatic and apply the scaling:
    double beta_internal = params_.compressibility / 0.0602214086;
    mu = 1.0 - (dt / params_.tau) * beta_internal * (params_.pressure - P_inst);

    if (mu <= 0.0) mu = 0.0;

    double scale = std::cbrt(mu);

    for (size_t i = 0; i < sys.n_atoms; ++i) {
        sys.pos_x[i] *= scale;
        sys.pos_y[i] *= scale;
        sys.pos_z[i] *= scale;
    }

    cell.scale(scale, scale, scale);
}
