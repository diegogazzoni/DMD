#include "force/force_engine.h"
#include "force/force_component.h"
#include "core/cell.h"
#include "core/system_data.h"
#include <gtest/gtest.h>

class MockComponent : public ForceComponent {
public:
    std::string name() const override { return "Mock"; }
    void compute(
        std::span<const double>, std::span<const double>, std::span<const double>,
        const Cell&,
        std::span<double>, std::span<double>, std::span<double>,
        double& energy_out) override
    {
        energy_out += 42.0;
    }
};

TEST(ForceEngineTest, SingleComponent) {
    ForceEngine engine;
    engine.add_component(std::make_unique<MockComponent>());

    SystemData sys(10);
    Cell cell(5.0, 5.0, 5.0);

    auto report = engine.compute(sys, cell);

    EXPECT_DOUBLE_EQ(report.total, 42.0);
    EXPECT_DOUBLE_EQ(sys.potential_energy, 42.0);
}

TEST(ForceEngineTest, TwoComponentsAccumulate) {
    ForceEngine engine;
    engine.add_component(std::make_unique<MockComponent>());
    engine.add_component(std::make_unique<MockComponent>());

    SystemData sys(10);
    Cell cell(5.0, 5.0, 5.0);

    auto report = engine.compute(sys, cell);

    EXPECT_DOUBLE_EQ(report.total, 84.0);
}

TEST(ForceEngineTest, ForcesAreZeroedBetweenCompute) {
    ForceEngine engine;
    engine.add_component(std::make_unique<MockComponent>());

    SystemData sys(5);
    sys.forces_x[2] = 100.0;

    Cell cell(5.0, 5.0, 5.0);
    engine.compute(sys, cell);

    EXPECT_DOUBLE_EQ(sys.forces_x[2], 100.0);
}
