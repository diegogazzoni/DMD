#pragma once

#include "thermostat.h"

struct NoseHooverParams {
    double temperature{300.0};
    double tau{0.1}; // relaxation time (ps)
};

class NoseHooverThermostat : public Thermostat {
public:
    explicit NoseHooverThermostat(NoseHooverParams params, size_t n_atoms);
    std::string name() const override;
    void apply(SystemData& sys, double dt) override;

    double zeta() const { return zeta_; }
    void set_zeta(double z) { zeta_ = z; }

private:
    NoseHooverParams params_;
    double zeta_{0.0};
    double Q_{1.0}; // Nose mass
};
