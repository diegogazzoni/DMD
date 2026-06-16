#pragma once

#include "thermostat.h"

struct BerendsenParams {
    double temperature{300.0};
    double tau{0.1}; // coupling time (ps)
};

class BerendsenThermostat : public Thermostat {
public:
    explicit BerendsenThermostat(BerendsenParams params);
    std::string name() const override;
    void apply(SystemData& sys, double dt) override;

private:
    BerendsenParams params_;
};
