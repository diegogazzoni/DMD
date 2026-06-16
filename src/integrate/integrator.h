#pragma once

#include <span>

struct SystemData;

class Integrator {
public:
    void step(
        SystemData& sys,
        double dt
    );
};
