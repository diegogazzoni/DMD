#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "sim/simulation_config.h"
#include "sim/simulation_engine.h"
#include "sim/json_config.h"
#include "sysbin/dmdin.h"
#include "core/system_data.h"
#include "core/cell.h"
#include <cmath>

namespace py = pybind11;

// Wrapper to expose SimulationEngine with a simpler interface
struct EngineWrapper {
    SimulationEngine engine_;
    double box_size_;
    size_t n_atoms_;

    EngineWrapper(SimulationConfig cfg)
        : engine_(build_simulation(std::move(cfg)))
        , box_size_(cfg.box_size)
        , n_atoms_(cfg.masses.size())
    {}

    void run() { engine_.run(); }

    void run_n(int steps) {
        for (int i = 0; i < steps; ++i)
            engine_.step();
    }

    double box_size() const { return box_size_; }
    size_t n_atoms() const { return n_atoms_; }

    py::array_t<double> forces() const {
        auto& sys = engine_.system();
        py::ssize_t n = static_cast<py::ssize_t>(sys.n_atoms);
        auto arr = py::array_t<double>({n, py::ssize_t(3)});
        auto buf = arr.template mutable_unchecked<2>();
        for (size_t i = 0; i < sys.n_atoms; ++i) {
            buf(i, 0) = sys.forces_x[i];
            buf(i, 1) = sys.forces_y[i];
            buf(i, 2) = sys.forces_z[i];
        }
        return arr;
    }

    int current_step() const { return engine_.current_step(); }
    double current_time() const { return engine_.current_time(); }

    // Access system data as numpy arrays
    py::array_t<double> positions() const {
        auto& sys = engine_.system();
        py::ssize_t n = static_cast<py::ssize_t>(sys.n_atoms);
        auto arr = py::array_t<double>({n, py::ssize_t(3)});
        auto buf = arr.template mutable_unchecked<2>();
        for (size_t i = 0; i < sys.n_atoms; ++i) {
            buf(i, 0) = sys.pos_x[i];
            buf(i, 1) = sys.pos_y[i];
            buf(i, 2) = sys.pos_z[i];
        }
        return arr;
    }

    py::array_t<double> velocities() const {
        auto& sys = engine_.system();
        py::ssize_t n = static_cast<py::ssize_t>(sys.n_atoms);
        auto arr = py::array_t<double>({n, py::ssize_t(3)});
        auto buf = arr.template mutable_unchecked<2>();
        for (size_t i = 0; i < sys.n_atoms; ++i) {
            buf(i, 0) = sys.vel_x[i];
            buf(i, 1) = sys.vel_y[i];
            buf(i, 2) = sys.vel_z[i];
        }
        return arr;
    }

    double potential_energy() const {
        return engine_.system().potential_energy;
    }
};

// Helper to build a SimulationConfig from Python-friendly structures
static SimulationConfig build_cfg_from_py(
    py::array_t<double> positions,
    py::array_t<double> masses,
    py::array_t<double> charges,
    py::array_t<int32_t> atom_types,
    double box_size,
    double dt,
    int n_steps,
    double lj_cutoff,
    double init_temperature,
    bool gen_vel,
    int seed)
{
    SimulationConfig cfg;
    cfg.box_size = box_size;
    cfg.dt = dt;
    cfg.n_steps = n_steps;
    cfg.lj_cutoff = lj_cutoff;
    cfg.use_lj = true;
    cfg.lj.n_types = 1;
    cfg.lj.sigma = {0.3405};
    cfg.lj.epsilon = {0.997};
    cfg.init_temperature = init_temperature;
    cfg.gen_vel = gen_vel;
    cfg.seed = seed;

    size_t n = static_cast<size_t>(positions.shape(0));
    auto pos = positions.unchecked<2>();
    for (size_t i = 0; i < n; ++i) {
        cfg.pos_x.push_back(pos(i, 0));
        cfg.pos_y.push_back(pos(i, 1));
        cfg.pos_z.push_back(pos(i, 2));
        cfg.vel_x.push_back(0.0);
        cfg.vel_y.push_back(0.0);
        cfg.vel_z.push_back(0.0);
    }

    auto m = masses.unchecked<1>();
    auto c = charges.unchecked<1>();
    auto t = atom_types.unchecked<1>();
    for (size_t i = 0; i < n; ++i) {
        cfg.masses.push_back(m(i));
        cfg.charges.push_back(c(i));
        cfg.atom_types.push_back(t(i));
    }

    return cfg;
}

PYBIND11_MODULE(_dmd_core, m) {
    m.doc() = "DMD molecular dynamics engine — C++ core bindings";

    // Helper to convert std::vector to py::array_t
    auto vec_to_arr_d = [](const std::vector<double>& v) {
        py::ssize_t sz = static_cast<py::ssize_t>(v.size());
        auto result = py::array_t<double>(sz, v.data());
        return result;
    };
    auto vec_to_arr_i = [](const std::vector<int>& v) {
        py::ssize_t sz = static_cast<py::ssize_t>(v.size());
        auto result = py::array_t<int>(sz, v.data());
        return result;
    };

    // ---- SimulationConfig binding ----
    py::class_<SimulationConfig>(m, "SimulationConfig")
        .def(py::init<>())
        .def_readwrite("dt", &SimulationConfig::dt)
        .def_readwrite("n_steps", &SimulationConfig::n_steps)
        .def_readwrite("box_size", &SimulationConfig::box_size)
        .def_readwrite("lj_cutoff", &SimulationConfig::lj_cutoff)
        .def_readwrite("use_lj", &SimulationConfig::use_lj)
        .def_readwrite("init_temperature", &SimulationConfig::init_temperature)
        .def_readwrite("gen_vel", &SimulationConfig::gen_vel)
        .def_readwrite("seed", &SimulationConfig::seed)
        .def_readwrite("trajectory_path", &SimulationConfig::trajectory_path)
        .def_readwrite("nstxout", &SimulationConfig::nstxout)
        .def_readwrite("nstvout", &SimulationConfig::nstvout)
        .def_readwrite("nstenergy", &SimulationConfig::nstenergy)
        .def_readwrite("energy_path", &SimulationConfig::energy_path)
        .def_readwrite("checkpoint_interval", &SimulationConfig::checkpoint_interval)
        .def_readwrite("checkpoint_path", &SimulationConfig::checkpoint_path)
        .def_readwrite("coulomb_type", &SimulationConfig::coulomb_type)
        .def_readwrite("coulomb_cutoff", &SimulationConfig::coulomb_cutoff)
        .def_readwrite("pme_order", &SimulationConfig::pme_order)
        .def_readwrite("pme_grid_spacing", &SimulationConfig::pme_grid_spacing)
        .def_readwrite("ewald_coeff", &SimulationConfig::ewald_coeff)
        .def_readwrite("constraint_type", &SimulationConfig::constraint_type)
        .def_readwrite("constraint_tolerance", &SimulationConfig::constraint_tolerance)
        // Vectors as numpy arrays
        .def_property("mass",        [&](SimulationConfig& c) { return vec_to_arr_d(c.masses); },
            [](SimulationConfig& c, py::array_t<double> a) {
                c.masses.assign(a.data(), a.data() + a.size());
            })
        .def_property("charges",     [&](SimulationConfig& c) { return vec_to_arr_d(c.charges); },
            [](SimulationConfig& c, py::array_t<double> a) {
                c.charges.assign(a.data(), a.data() + a.size());
            })
        .def_property("pos_x",       [&](SimulationConfig& c) { return vec_to_arr_d(c.pos_x); },
            [](SimulationConfig& c, py::array_t<double> a) {
                c.pos_x.assign(a.data(), a.data() + a.size());
            })
        .def_property("pos_y",       [&](SimulationConfig& c) { return vec_to_arr_d(c.pos_y); },
            [](SimulationConfig& c, py::array_t<double> a) {
                c.pos_y.assign(a.data(), a.data() + a.size());
            })
        .def_property("pos_z",       [&](SimulationConfig& c) { return vec_to_arr_d(c.pos_z); },
            [](SimulationConfig& c, py::array_t<double> a) {
                c.pos_z.assign(a.data(), a.data() + a.size());
            })
        .def_property("vel_x",       [&](SimulationConfig& c) { return vec_to_arr_d(c.vel_x); },
            [](SimulationConfig& c, py::array_t<double> a) {
                c.vel_x.assign(a.data(), a.data() + a.size());
            })
        .def_property("vel_y",       [&](SimulationConfig& c) { return vec_to_arr_d(c.vel_y); },
            [](SimulationConfig& c, py::array_t<double> a) {
                c.vel_y.assign(a.data(), a.data() + a.size());
            })
        .def_property("vel_z",       [&](SimulationConfig& c) { return vec_to_arr_d(c.vel_z); },
            [](SimulationConfig& c, py::array_t<double> a) {
                c.vel_z.assign(a.data(), a.data() + a.size());
            })
        .def_property("atom_types",  [&](SimulationConfig& c) { return vec_to_arr_i(c.atom_types); },
            [](SimulationConfig& c, py::array_t<int> a) {
                c.atom_types.assign(a.data(), a.data() + a.size());
            })
        // LJ params
        .def_property("lj_sigma",    [&](SimulationConfig& c) { return vec_to_arr_d(c.lj.sigma); },
            [](SimulationConfig& c, py::array_t<double> a) {
                c.lj.sigma.assign(a.data(), a.data() + a.size());
                if (!c.lj.sigma.empty())
                    c.lj.n_types = static_cast<size_t>(std::sqrt(c.lj.sigma.size()));
            })
        .def_property("lj_epsilon",  [&](SimulationConfig& c) { return vec_to_arr_d(c.lj.epsilon); },
            [](SimulationConfig& c, py::array_t<double> a) {
                c.lj.epsilon.assign(a.data(), a.data() + a.size());
                if (!c.lj.epsilon.empty())
                    c.lj.n_types = static_cast<size_t>(std::sqrt(c.lj.epsilon.size()));
            })
        // Bonds (empty if no bonds)
        .def_property("bond_i",
            [&](SimulationConfig& c) {
                if (c.bonds.empty()) return py::array_t<int>();
                return vec_to_arr_i(c.bonds[0].i);
            },
            [](SimulationConfig& c, py::array_t<int> a) {
                if (c.bonds.empty()) c.bonds.resize(1);
                c.bonds[0].i.assign(a.data(), a.data() + a.size());
                c.bonds[0].n = a.size();
            })
        .def_property("bond_j",
            [&](SimulationConfig& c) {
                if (c.bonds.empty()) return py::array_t<int>();
                return vec_to_arr_i(c.bonds[0].j);
            },
            [](SimulationConfig& c, py::array_t<int> a) {
                if (c.bonds.empty()) c.bonds.resize(1);
                c.bonds[0].j.assign(a.data(), a.data() + a.size());
            })
        .def_property("bond_k",
            [&](SimulationConfig& c) {
                if (c.bonds.empty()) return py::array_t<double>();
                return vec_to_arr_d(c.bonds[0].k);
            },
            [](SimulationConfig& c, py::array_t<double> a) {
                if (c.bonds.empty()) c.bonds.resize(1);
                c.bonds[0].k.assign(a.data(), a.data() + a.size());
            })
        .def_property("bond_r0",
            [&](SimulationConfig& c) {
                if (c.bonds.empty()) return py::array_t<double>();
                return vec_to_arr_d(c.bonds[0].r0);
            },
            [](SimulationConfig& c, py::array_t<double> a) {
                if (c.bonds.empty()) c.bonds.resize(1);
                c.bonds[0].r0.assign(a.data(), a.data() + a.size());
            });

    // ---- Engine wrapper ----
    py::class_<EngineWrapper>(m, "Engine")
        .def(py::init<SimulationConfig>(), py::arg("cfg"))
        .def("run", &EngineWrapper::run)
        .def("run_n", &EngineWrapper::run_n, py::arg("steps"))
        .def_property_readonly("current_step", &EngineWrapper::current_step)
        .def_property_readonly("current_time", &EngineWrapper::current_time)
        .def_property_readonly("positions", &EngineWrapper::positions)
        .def_property_readonly("velocities", &EngineWrapper::velocities)
        .def_property_readonly("forces", &EngineWrapper::forces)
        .def_property_readonly("box_size", &EngineWrapper::box_size)
        .def_property_readonly("n_atoms", &EngineWrapper::n_atoms)
        .def_property_readonly("potential_energy", &EngineWrapper::potential_energy);

    // ---- Free functions ----
    m.def("read_dmdin", &read_dmdin, "Read a .dmdin binary file",
          py::arg("path"));

    m.def("write_dmdin", &write_dmdin, "Write a .dmdin binary file",
          py::arg("path"), py::arg("cfg"));

    m.def("apply_json_config", &apply_json_config,
          "Apply a config.json to a SimulationConfig",
          py::arg("cfg"), py::arg("path"));

    m.def("generate_config_template", &generate_config_template,
          "Generate a config.json template string");

    m.def("build_cfg",
          &build_cfg_from_py,
          "Build a SimulationConfig from numpy arrays",
          py::arg("positions"), py::arg("masses"),
          py::arg("charges"), py::arg("atom_types"),
          py::arg("box_size") = 3.0,
          py::arg("dt") = 0.002,
          py::arg("n_steps") = 1000,
          py::arg("lj_cutoff") = 1.2,
          py::arg("init_temperature") = 300.0,
          py::arg("gen_vel") = false,
          py::arg("seed") = 42);
}
