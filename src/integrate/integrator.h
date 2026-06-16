#pragma once

struct SystemData;

class Integrator {
public:
    void half_kick(SystemData& sys, double dt);
    void advance(SystemData& sys, double dt);
};
