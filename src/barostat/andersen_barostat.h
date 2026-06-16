#pragma once

#include "barostat.h"

struct AndersenBarostatParams {
    double pressure{1.0};  // target pressure (bar)
    double tau{1.0};       // relaxation time (ps)
};

class AndersenBarostat : public Barostat {
public:
    AndersenBarostat(AndersenBarostatParams params, size_t n_atoms);
    std::string name() const override;
    void apply(SystemData& sys, Cell& cell, double dt) override;

    double piston_velocity() const { return eps_dot_; }
    void set_piston_velocity(double v) { eps_dot_ = v; }

private:
    AndersenBarostatParams params_;
    double eps_dot_{0.0}; // strain rate (piston velocity)
    double W_{1.0};       // piston mass
};
