#include "checkpoint.h"
#include "core/system_data.h"
#include <fstream>
#include <cstdint>

static constexpr uint32_t MAGIC = 0x444D4450;
static constexpr uint32_t VERSION = 1;

void checkpoint::save(const std::string& path, const SystemData& sys, int lj_step_since_rebuild) {
    std::ofstream os(path, std::ios::binary);
    uint64_t n = sys.n_atoms;

    os.write(reinterpret_cast<const char*>(&MAGIC), sizeof(MAGIC));
    os.write(reinterpret_cast<const char*>(&VERSION), sizeof(VERSION));
    os.write(reinterpret_cast<const char*>(&n), sizeof(n));
    os.write(reinterpret_cast<const char*>(&sys.step), sizeof(sys.step));
    os.write(reinterpret_cast<const char*>(&sys.time), sizeof(sys.time));
    os.write(reinterpret_cast<const char*>(&sys.potential_energy), sizeof(sys.potential_energy));

    auto write_vec = [&](const auto& v) {
        os.write(reinterpret_cast<const char*>(v.data()), n * sizeof(v[0]));
    };
    write_vec(sys.pos_x);
    write_vec(sys.pos_y);
    write_vec(sys.pos_z);
    write_vec(sys.vel_x);
    write_vec(sys.vel_y);
    write_vec(sys.vel_z);
    write_vec(sys.forces_x);
    write_vec(sys.forces_y);
    write_vec(sys.forces_z);
    write_vec(sys.masses);
    write_vec(sys.charges);

    os.write(reinterpret_cast<const char*>(sys.atom_types.data()), n * sizeof(int));

    os.write(reinterpret_cast<const char*>(&lj_step_since_rebuild), sizeof(lj_step_since_rebuild));
}

void checkpoint::load(const std::string& path, SystemData& sys, int& lj_step_since_rebuild) {
    std::ifstream is(path, std::ios::binary);
    uint32_t magic, version;
    uint64_t n;

    is.read(reinterpret_cast<char*>(&magic), sizeof(magic));
    is.read(reinterpret_cast<char*>(&version), sizeof(version));
    is.read(reinterpret_cast<char*>(&n), sizeof(n));

    if (n != sys.n_atoms) {
        sys = SystemData(static_cast<size_t>(n));
    }

    is.read(reinterpret_cast<char*>(&sys.step), sizeof(sys.step));
    is.read(reinterpret_cast<char*>(&sys.time), sizeof(sys.time));
    is.read(reinterpret_cast<char*>(&sys.potential_energy), sizeof(sys.potential_energy));

    auto read_vec = [&](auto& v) {
        is.read(reinterpret_cast<char*>(v.data()), n * sizeof(v[0]));
    };
    read_vec(sys.pos_x);
    read_vec(sys.pos_y);
    read_vec(sys.pos_z);
    read_vec(sys.vel_x);
    read_vec(sys.vel_y);
    read_vec(sys.vel_z);
    read_vec(sys.forces_x);
    read_vec(sys.forces_y);
    read_vec(sys.forces_z);
    read_vec(sys.masses);
    read_vec(sys.charges);

    is.read(reinterpret_cast<char*>(sys.atom_types.data()), n * sizeof(int));

    is.read(reinterpret_cast<char*>(&lj_step_since_rebuild), sizeof(lj_step_since_rebuild));
}
