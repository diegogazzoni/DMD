#pragma once

#include <vector>

struct Cell;
struct SystemData;

struct ConstraintsParams {
    std::vector<int> i, j;
    std::vector<double> distance_target;
    double tolerance{1e-8};
};

class Constraints {
public:
    explicit Constraints(ConstraintsParams params);
    void apply(SystemData& sys, const Cell& cell);
    void apply_rattle(SystemData& sys, const Cell& cell);

private:
    ConstraintsParams params_;
};
