#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "sim/simulation_config.h"
#include "sim/simulation_engine.h"
#include "core/system_data.h"
#include "core/cell.h"
#include "trajectory/h5md_writer.h"
#include <cmath>

namespace py = pybind11;

// Wrapper to expose SimulationEngine with a simpler interface
struct EngineWrapper {
    SimulationEngine engine_;
    std::unique_ptr<H5MDWriter> trajectory_writer_;
    double box_size_;
    size_t n_atoms_;

    EngineWrapper(SimulationConfig cfg)
        : engine_(build_simulation(std::move(cfg)))
        , box_size_(cfg.box_size)
        , n_atoms_(cfg.masses.size())
    {}

    void run() { engine_.run(); }

    void run_n(int steps) {
        if (engine_.config().trajectory_interval > 0 && !trajectory_writer_) {
            auto& sys = engine_.system();
            trajectory_writer_ = std::make_unique<H5MDWriter>(
                engine_.config().trajectory_path, sys.n_atoms,
                /*write_vel=*/true);
        }
        for (int i = 0; i < steps; ++i) {
            engine_.step();
            if (trajectory_writer_ &&
                engine_.current_step() % engine_.config().trajectory_interval == 0) {
                auto& sys = engine_.system();
                trajectory_writer_->write_frame(sys, engine_.cell());
            }
        }
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

    // ---- Nested config types ----
    py::class_<ThermostatConfig>(m, "ThermostatConfig")
        .def(py::init<>())
        .def_readwrite("type", &ThermostatConfig::type)
        .def_readwrite("temperature", &ThermostatConfig::temperature)
        .def_readwrite("tau", &ThermostatConfig::tau)
        .def_readwrite("frequency", &ThermostatConfig::frequency);

    py::class_<BarostatConfig>(m, "BarostatConfig")
        .def(py::init<>())
        .def_readwrite("type", &BarostatConfig::type)
        .def_readwrite("pressure", &BarostatConfig::pressure)
        .def_readwrite("tau", &BarostatConfig::tau)
        .def_readwrite("compressibility", &BarostatConfig::compressibility);

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
        .def_readwrite("constraint_i", &SimulationConfig::constraint_i)
        .def_readwrite("constraint_j", &SimulationConfig::constraint_j)
        .def_readwrite("constraint_dist", &SimulationConfig::constraint_dist)
        .def_readwrite("thermostat", &SimulationConfig::thermostat)
        .def_readwrite("barostat", &SimulationConfig::barostat)
        .def_readwrite("excl_i", &SimulationConfig::excl_i)
        .def_readwrite("excl_j", &SimulationConfig::excl_j)
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
            })
        // Angles (empty if no angles)
        .def_property("angle_i",
            [&](SimulationConfig& c) {
                if (c.angles.empty()) return py::array_t<int>();
                return vec_to_arr_i(c.angles[0].i);
            },
            [](SimulationConfig& c, py::array_t<int> a) {
                if (c.angles.empty()) c.angles.resize(1);
                c.angles[0].i.assign(a.data(), a.data() + a.size());
                c.angles[0].n = a.size();
            })
        .def_property("angle_j",
            [&](SimulationConfig& c) {
                if (c.angles.empty()) return py::array_t<int>();
                return vec_to_arr_i(c.angles[0].j);
            },
            [](SimulationConfig& c, py::array_t<int> a) {
                if (c.angles.empty()) c.angles.resize(1);
                c.angles[0].j.assign(a.data(), a.data() + a.size());
            })
        .def_property("angle_k",
            [&](SimulationConfig& c) {
                if (c.angles.empty()) return py::array_t<int>();
                return vec_to_arr_i(c.angles[0].k);
            },
            [](SimulationConfig& c, py::array_t<int> a) {
                if (c.angles.empty()) c.angles.resize(1);
                c.angles[0].k.assign(a.data(), a.data() + a.size());
            })
        .def_property("angle_k_theta",
            [&](SimulationConfig& c) {
                if (c.angles.empty()) return py::array_t<double>();
                return vec_to_arr_d(c.angles[0].k_theta);
            },
            [](SimulationConfig& c, py::array_t<double> a) {
                if (c.angles.empty()) c.angles.resize(1);
                c.angles[0].k_theta.assign(a.data(), a.data() + a.size());
            })
        .def_property("angle_theta0",
            [&](SimulationConfig& c) {
                if (c.angles.empty()) return py::array_t<double>();
                return vec_to_arr_d(c.angles[0].theta0);
            },
            [](SimulationConfig& c, py::array_t<double> a) {
                if (c.angles.empty()) c.angles.resize(1);
                c.angles[0].theta0.assign(a.data(), a.data() + a.size());
            })
        // Dihedrals (empty if no dihedrals)
        .def_property("dihedral_i",
            [&](SimulationConfig& c) {
                if (c.dihedrals.empty()) return py::array_t<int>();
                return vec_to_arr_i(c.dihedrals[0].i);
            },
            [](SimulationConfig& c, py::array_t<int> a) {
                if (c.dihedrals.empty()) c.dihedrals.resize(1);
                c.dihedrals[0].i.assign(a.data(), a.data() + a.size());
                c.dihedrals[0].n = a.size();
            })
        .def_property("dihedral_j",
            [&](SimulationConfig& c) {
                if (c.dihedrals.empty()) return py::array_t<int>();
                return vec_to_arr_i(c.dihedrals[0].j);
            },
            [](SimulationConfig& c, py::array_t<int> a) {
                if (c.dihedrals.empty()) c.dihedrals.resize(1);
                c.dihedrals[0].j.assign(a.data(), a.data() + a.size());
            })
        .def_property("dihedral_k",
            [&](SimulationConfig& c) {
                if (c.dihedrals.empty()) return py::array_t<int>();
                return vec_to_arr_i(c.dihedrals[0].k);
            },
            [](SimulationConfig& c, py::array_t<int> a) {
                if (c.dihedrals.empty()) c.dihedrals.resize(1);
                c.dihedrals[0].k.assign(a.data(), a.data() + a.size());
            })
        .def_property("dihedral_l",
            [&](SimulationConfig& c) {
                if (c.dihedrals.empty()) return py::array_t<int>();
                return vec_to_arr_i(c.dihedrals[0].l);
            },
            [](SimulationConfig& c, py::array_t<int> a) {
                if (c.dihedrals.empty()) c.dihedrals.resize(1);
                c.dihedrals[0].l.assign(a.data(), a.data() + a.size());
            })
        .def_property("dihedral_k_phi",
            [&](SimulationConfig& c) {
                if (c.dihedrals.empty()) return py::array_t<double>();
                return vec_to_arr_d(c.dihedrals[0].k_phi);
            },
            [](SimulationConfig& c, py::array_t<double> a) {
                if (c.dihedrals.empty()) c.dihedrals.resize(1);
                c.dihedrals[0].k_phi.assign(a.data(), a.data() + a.size());
            })
        .def_property("dihedral_periodicity",
            [&](SimulationConfig& c) {
                if (c.dihedrals.empty()) return py::array_t<int>();
                return vec_to_arr_i(c.dihedrals[0].periodicity);
            },
            [](SimulationConfig& c, py::array_t<int> a) {
                if (c.dihedrals.empty()) c.dihedrals.resize(1);
                c.dihedrals[0].periodicity.assign(a.data(), a.data() + a.size());
            })
        .def_property("dihedral_phi0",
            [&](SimulationConfig& c) {
                if (c.dihedrals.empty()) return py::array_t<double>();
                return vec_to_arr_d(c.dihedrals[0].phi0);
            },
            [](SimulationConfig& c, py::array_t<double> a) {
                if (c.dihedrals.empty()) c.dihedrals.resize(1);
                c.dihedrals[0].phi0.assign(a.data(), a.data() + a.size());
            })
        // Impropers (empty if no impropers)
        .def_property("improper_i",
            [&](SimulationConfig& c) {
                if (c.impropers.empty()) return py::array_t<int>();
                return vec_to_arr_i(c.impropers[0].i);
            },
            [](SimulationConfig& c, py::array_t<int> a) {
                if (c.impropers.empty()) c.impropers.resize(1);
                c.impropers[0].i.assign(a.data(), a.data() + a.size());
                c.impropers[0].n = a.size();
            })
        .def_property("improper_j",
            [&](SimulationConfig& c) {
                if (c.impropers.empty()) return py::array_t<int>();
                return vec_to_arr_i(c.impropers[0].j);
            },
            [](SimulationConfig& c, py::array_t<int> a) {
                if (c.impropers.empty()) c.impropers.resize(1);
                c.impropers[0].j.assign(a.data(), a.data() + a.size());
            })
        .def_property("improper_k",
            [&](SimulationConfig& c) {
                if (c.impropers.empty()) return py::array_t<int>();
                return vec_to_arr_i(c.impropers[0].k);
            },
            [](SimulationConfig& c, py::array_t<int> a) {
                if (c.impropers.empty()) c.impropers.resize(1);
                c.impropers[0].k.assign(a.data(), a.data() + a.size());
            })
        .def_property("improper_l",
            [&](SimulationConfig& c) {
                if (c.impropers.empty()) return py::array_t<int>();
                return vec_to_arr_i(c.impropers[0].l);
            },
            [](SimulationConfig& c, py::array_t<int> a) {
                if (c.impropers.empty()) c.impropers.resize(1);
                c.impropers[0].l.assign(a.data(), a.data() + a.size());
            })
        .def_property("improper_k_phi",
            [&](SimulationConfig& c) {
                if (c.impropers.empty()) return py::array_t<double>();
                return vec_to_arr_d(c.impropers[0].k_phi);
            },
            [](SimulationConfig& c, py::array_t<double> a) {
                if (c.impropers.empty()) c.impropers.resize(1);
                c.impropers[0].k_phi.assign(a.data(), a.data() + a.size());
            })
        .def_property("improper_phi0",
            [&](SimulationConfig& c) {
                if (c.impropers.empty()) return py::array_t<double>();
                return vec_to_arr_d(c.impropers[0].phi0);
            },
            [](SimulationConfig& c, py::array_t<double> a) {
                if (c.impropers.empty()) c.impropers.resize(1);
                c.impropers[0].phi0.assign(a.data(), a.data() + a.size());
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

}
