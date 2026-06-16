#pragma once

#include <span>
#include <vector>

struct SystemData;

enum class ConstraintAlgorithm {
    shake,
    rattle
};

struct ConstraintsParams {
    size_t n{0};
    std::vector<int> i, j;
    std::vector<double> distance_target;
    ConstraintAlgorithm algorithm{ConstraintAlgorithm::shake};
};

class Constraints {
public:
    explicit Constraints(ConstraintsParams params);
    void apply(SystemData& sys, double dt);

private:
    ConstraintsParams params_;
};
