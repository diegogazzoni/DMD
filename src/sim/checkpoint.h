#pragma once

#include <string>

struct SystemData;
class ForceEngine;
class LennardJones;

namespace checkpoint {

void save(const std::string& path, const SystemData& sys, int lj_step_since_rebuild);

void load(const std::string& path, SystemData& sys, int& lj_step_since_rebuild);

} // namespace checkpoint
