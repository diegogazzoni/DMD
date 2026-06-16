#pragma once

#include "barostat.h"

struct BerendsenBarostatParams {
    double pressure{1.0};       // target pressure (bar)
    double tau{1.0};            // coupling time (ps)
    double compressibility{4.5e-5}; // isothermal compressibility (bar^-1)
};

class BerendsenBarostat : public Barostat {
public:
    explicit BerendsenBarostat(BerendsenBarostatParams params);
    std::string name() const override;
    void apply(SystemData& sys, Cell& cell, double dt) override;

private:
    BerendsenBarostatParams params_;
};
