#pragma once

#include "simulation_config.h"
#include <string>

void apply_json_config(SimulationConfig& cfg, const std::string& json_path);
std::string generate_config_template();
