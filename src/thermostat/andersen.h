#pragma once

#include "thermostat.h"

struct AndersenParams {
    double temperature{300.0};
    double frequency{10.0}; // collision frequency (ps^-1)
};

class AndersenThermostat : public Thermostat {
public:
    explicit AndersenThermostat(AndersenParams params);
    std::string name() const override;
    void apply(SystemData& sys, double dt) override;

private:
    AndersenParams params_;
    unsigned int seed_{42};
};
