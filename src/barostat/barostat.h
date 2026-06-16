#pragma once

#include <string>

struct SystemData;
struct Cell;

double virial(const SystemData& sys);
double pressure(double ekin, double virial_val, double volume);
double pressure(const SystemData& sys, double volume);

class Barostat {
public:
    virtual ~Barostat() = default;
    virtual std::string name() const = 0;
    virtual void apply(SystemData& sys, Cell& cell, double dt) = 0;
};
