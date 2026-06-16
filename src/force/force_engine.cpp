#include "force_engine.h"
#include "force/force_component.h"
#include "core/cell.h"
#include "core/system_data.h"

void ForceEngine::add_component(std::unique_ptr<ForceComponent> component) {
    components_.push_back(std::move(component));
}

EnergyReport ForceEngine::compute(SystemData& sys, const Cell& cell) {
    EnergyReport report;
    report.total = 0.0;
    sys.potential_energy = 0.0;

    for (auto& comp : components_) {
        comp->compute(
            sys.x(), sys.y(), sys.z(),
            cell,
            sys.fx(), sys.fy(), sys.fz(),
            sys.potential_energy
        );
    }
    report.total = sys.potential_energy;
    return report;
}
