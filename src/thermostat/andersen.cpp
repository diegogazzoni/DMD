#include "andersen.h"
#include "core/system_data.h"
#include "core/constants.h"
#include <cmath>
#include <random>

AndersenThermostat::AndersenThermostat(AndersenParams params)
    : params_(std::move(params)) {}

std::string AndersenThermostat::name() const {
    return "Andersen";
}

void AndersenThermostat::apply(SystemData& sys, double dt) {
    double prob = params_.frequency * dt;
    if (prob <= 0.0) return;

    double sigma_v = std::sqrt(dmd_constants::KB * params_.temperature / 39.948); // placeholder mass for scaling
    std::mt19937_64 rng(seed_);
    std::uniform_real_distribution<double> uni(0.0, 1.0);
    std::normal_distribution<double> ndist(0.0, 1.0);

    for (size_t i = 0; i < sys.n_atoms; ++i) {
        if (uni(rng) < prob) {
            double sigma = std::sqrt(dmd_constants::KB * params_.temperature / sys.masses[i]);
            sys.vel_x[i] = ndist(rng) * sigma;
            sys.vel_y[i] = ndist(rng) * sigma;
            sys.vel_z[i] = ndist(rng) * sigma;
        }
    }
    ++seed_;
}
