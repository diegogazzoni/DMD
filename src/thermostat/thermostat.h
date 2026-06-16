#pragma once

#include <string>

struct SystemData;

double kinetic_energy(const SystemData& sys);
double temperature(double ekin, size_t n_atoms);
double temperature(const SystemData& sys);

class Thermostat {
public:
    virtual ~Thermostat() = default;
    virtual std::string name() const = 0;
    virtual void apply(SystemData& sys, double dt) = 0;
};
