#include "json_config.h"
#include "simulation_config.h"
#include <nlohmann/json.hpp>
#include <fstream>
#include <stdexcept>
#include <sstream>
#include <set>
#include <string>

using json = nlohmann::json;

// ---- Schema ----

// The canonical schema: all keys, all values are type-placeholders.
// Every key in user JSON must appear here; every key here must appear in user JSON.
static const json& schema() {
    static const json s = json::object({
        {"run", json::object({
            {"dt", 0.0},
            {"n_steps", 0},
            {"init_temperature", 0.0},
            {"gen_vel", false},
            {"seed", 0},
        })},
        {"output", json::object({
            {"trajectory_path", ""},
            {"nstxout", 0},
            {"nstvout", 0},
            {"nstenergy", 0},
            {"energy_path", ""},
            {"checkpoint_interval", 0},
            {"checkpoint_path", ""},
        })},
        {"lj", json::object({
            {"cutoff", 0.0},
        })},
        {"electrostatics", json::object({
            {"coulomb_type", ""},
            {"cutoff", 0.0},
            {"pme_order", 0},
            {"pme_grid_spacing", 0.0},
            {"ewald_coeff", 0.0},
        })},
        {"thermostat", json::object({
            {"type", ""},
            {"temperature", 0.0},
            {"tau", 0.0},
            {"frequency", 0.0},
        })},
        {"barostat", json::object({
            {"type", ""},
            {"pressure", 0.0},
            {"tau", 0.0},
            {"compressibility", 0.0},
        })},
        {"constraints", json::object({
            {"type", ""},
            {"tolerance", 0.0},
        })},
    });
    return s;
}

// ---- Validation ----

static void collect_keys(const json& j, const std::string& prefix,
                         std::set<std::string>& out) {
    if (j.is_object()) {
        for (auto& [k, v] : j.items()) {
            std::string path = prefix.empty() ? k : prefix + "." + k;
            out.insert(path);
            collect_keys(v, path, out);
        }
    }
}

static std::string json_type_name(const json& j) {
    if (j.is_null()) return "null";
    if (j.is_boolean()) return "boolean";
    if (j.is_number_integer()) return "integer";
    if (j.is_number_float()) return "float";
    if (j.is_string()) return "string";
    if (j.is_array()) return "array";
    if (j.is_object()) return "object";
    return "unknown";
}

static bool json_type_compatible(const json& schema_val, const json& input_val) {
    if (schema_val.is_number_float() && input_val.is_number()) return true;
    if (schema_val.is_number_integer() && input_val.is_number_integer()) return true;
    return schema_val.type() == input_val.type();
}

static void validate_schema(const json& input, const json& sch,
                            const std::string& path) {
    // Check unknown keys
    for (auto& [k, v] : input.items()) {
        std::string child_path = path.empty() ? k : path + "." + k;
        if (!sch.contains(k)) {
            throw std::runtime_error("Unknown key '" + child_path + "'");
        }
        if (v.is_object() && sch[k].is_object()) {
            validate_schema(v, sch[k], child_path);
        } else if (json_type_compatible(sch[k], v)) {
            // type OK
        } else {
            throw std::runtime_error("Type mismatch for '" + child_path +
                "': expected " + json_type_name(sch[k]) +
                ", got " + json_type_name(v));
        }
    }
    // Check missing keys
    for (auto& [k, v] : sch.items()) {
        std::string child_path = path.empty() ? k : path + "." + k;
        if (!input.contains(k)) {
            throw std::runtime_error("Missing required key '" + child_path + "'");
        }
        if (v.is_object() && input[k].is_object()) {
            validate_schema(input[k], v, child_path);
        }
    }
}

// ---- Reader ----

static double read_double(const json& j, const std::string& section,
                           const std::string& key) {
    try { return j[key].get<double>(); }
    catch (...) { throw std::runtime_error("Cannot read '" + section + "." + key + "' as number"); }
}
static int read_int(const json& j, const std::string& section,
                     const std::string& key) {
    try { return j[key].get<int>(); }
    catch (...) { throw std::runtime_error("Cannot read '" + section + "." + key + "' as integer"); }
}
static std::string read_string(const json& j, const std::string& section,
                                const std::string& key) {
    try { return j[key].get<std::string>(); }
    catch (...) { throw std::runtime_error("Cannot read '" + section + "." + key + "' as string"); }
}
static bool read_bool(const json& j, const std::string& section,
                       const std::string& key) {
    try { return j[key].get<bool>(); }
    catch (...) { throw std::runtime_error("Cannot read '" + section + "." + key + "' as boolean"); }
}

void apply_json_config(SimulationConfig& cfg, const std::string& json_path) {
    std::ifstream is(json_path);
    if (!is) throw std::runtime_error("Cannot open config file: " + json_path);

    json j;
    try { is >> j; }
    catch (const json::parse_error& e) {
        throw std::runtime_error("JSON parse error in " + json_path + ": " + e.what());
    }

    // Validate against schema
    const json& sch = schema();
    validate_schema(j, sch, "");

    // ---- run ----
    auto& run = j["run"];
    cfg.dt = read_double(run, "run", "dt");
    cfg.n_steps = read_int(run, "run", "n_steps");
    cfg.init_temperature = read_double(run, "run", "init_temperature");
    cfg.gen_vel = read_bool(run, "run", "gen_vel");
    cfg.seed = read_int(run, "run", "seed");

    // ---- output ----
    auto& out = j["output"];
    cfg.trajectory_path = read_string(out, "output", "trajectory_path");
    cfg.nstxout = read_int(out, "output", "nstxout");
    cfg.nstvout = read_int(out, "output", "nstvout");
    cfg.nstenergy = read_int(out, "output", "nstenergy");
    cfg.energy_path = read_string(out, "output", "energy_path");
    cfg.checkpoint_interval = read_int(out, "output", "checkpoint_interval");
    cfg.checkpoint_path = read_string(out, "output", "checkpoint_path");

    // ---- lj ----
    auto& lj = j["lj"];
    cfg.lj_cutoff = read_double(lj, "lj", "cutoff");
    cfg.use_lj = true;

    // ---- electrostatics ----
    auto& es = j["electrostatics"];
    cfg.coulomb_type = read_string(es, "electrostatics", "coulomb_type");
    cfg.coulomb_cutoff = read_double(es, "electrostatics", "cutoff");
    cfg.pme_order = read_int(es, "electrostatics", "pme_order");
    cfg.pme_grid_spacing = read_double(es, "electrostatics", "pme_grid_spacing");
    cfg.ewald_coeff = read_double(es, "electrostatics", "ewald_coeff");

    // ---- thermostat ----
    auto& tc = j["thermostat"];
    cfg.thermostat.type = read_string(tc, "thermostat", "type");
    cfg.thermostat.temperature = read_double(tc, "thermostat", "temperature");
    cfg.thermostat.tau = read_double(tc, "thermostat", "tau");
    cfg.thermostat.frequency = read_double(tc, "thermostat", "frequency");

    // ---- barostat ----
    auto& bc = j["barostat"];
    cfg.barostat.type = read_string(bc, "barostat", "type");
    cfg.barostat.pressure = read_double(bc, "barostat", "pressure");
    cfg.barostat.tau = read_double(bc, "barostat", "tau");
    cfg.barostat.compressibility = read_double(bc, "barostat", "compressibility");

    // ---- constraints ----
    auto& cc = j["constraints"];
    cfg.constraint_type = read_string(cc, "constraints", "type");
    cfg.constraint_tolerance = read_double(cc, "constraints", "tolerance");
}

// ---- Template ----

std::string generate_config_template() {
    // Build a pretty-printed JSON with example values
    json t;
    t["run"]["dt"] = 0.002;
    t["run"]["n_steps"] = 10000;
    t["run"]["init_temperature"] = 300.0;
    t["run"]["gen_vel"] = false;
    t["run"]["seed"] = 42;

    t["output"]["trajectory_path"] = "trajectory.h5";
    t["output"]["nstxout"] = 100;
    t["output"]["nstvout"] = 100;
    t["output"]["nstenergy"] = 10;
    t["output"]["energy_path"] = "energy.h5";
    t["output"]["checkpoint_interval"] = 500;
    t["output"]["checkpoint_path"] = "checkpoint.bin";

    t["lj"]["cutoff"] = 1.2;

    t["electrostatics"]["coulomb_type"] = "cutoff";
    t["electrostatics"]["cutoff"] = 1.2;
    t["electrostatics"]["pme_order"] = 4;
    t["electrostatics"]["pme_grid_spacing"] = 0.12;
    t["electrostatics"]["ewald_coeff"] = 0.34;

    t["thermostat"]["type"] = "none";
    t["thermostat"]["temperature"] = 300.0;
    t["thermostat"]["tau"] = 0.1;
    t["thermostat"]["frequency"] = 10.0;

    t["barostat"]["type"] = "none";
    t["barostat"]["pressure"] = 1.0;
    t["barostat"]["tau"] = 1.0;
    t["barostat"]["compressibility"] = 4.5e-5;

    t["constraints"]["type"] = "none";
    t["constraints"]["tolerance"] = 1e-6;

    return t.dump(2);
}
