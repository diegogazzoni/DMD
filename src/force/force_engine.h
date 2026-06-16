#pragma once

#include <memory>
#include <vector>

class ForceComponent;
struct Cell;
struct SystemData;

struct EnergyReport {
    double total{0.0};
};

class ForceEngine {
public:
    void add_component(std::unique_ptr<ForceComponent> component);

    EnergyReport compute(SystemData& sys, const Cell& cell);

    const std::vector<std::unique_ptr<ForceComponent>>& components() const {
        return components_;
    }

private:
    std::vector<std::unique_ptr<ForceComponent>> components_;
};
