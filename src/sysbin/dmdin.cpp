#include "dmdin.h"
#include "core/cell.h"
#include <fstream>
#include <vector>
#include <cstdint>
#include <cmath>

static constexpr uint32_t DMDIN_MAGIC = 0x444D444E;
static constexpr uint32_t DMDIN_VERSION = 1;
static constexpr uint32_t POS_HAS_VEL = 1;

struct DMDinHeader {
    uint32_t magic;
    uint32_t version;
    uint64_t n_atoms;
    uint32_t cell_type;
    int32_t  n_types;
    uint32_t pos_flags;
    uint32_t reserved;
    uint64_t offset_ff;
    uint64_t offset_pos;
    double   box[9];
};

static_assert(sizeof(DMDinHeader) == 120, "DMDinHeader must be 120 bytes");

static int cell_type_to_enum(const char* name) {
    std::string s(name);
    if (s == "orthorhombic") return 0;
    if (s == "triclinic") return 1;
    if (s == "hexagonal") return 2;
    if (s == "rhombic_dodecahedron") return 3;
    return 0;
}

static CellType int_to_cell_type(uint32_t t) {
    switch (t) {
        case 1: return CellType::triclinic;
        case 2: return CellType::hexagonal;
        case 3: return CellType::rhombic_dodecahedron;
        default: return CellType::orthorhombic;
    }
}

// ---- Reader ----

SimulationConfig read_dmdin(const std::string& path) {
    std::ifstream is(path, std::ios::binary);
    if (!is) throw std::runtime_error("Cannot open " + path);

    DMDinHeader header;
    is.read(reinterpret_cast<char*>(&header), sizeof(header));
    if (header.magic != DMDIN_MAGIC)
        throw std::runtime_error("Bad magic in " + path);
    if (header.version != DMDIN_VERSION)
        throw std::runtime_error("Unsupported version " + std::to_string(header.version));

    SimulationConfig cfg;
    size_t n = static_cast<size_t>(header.n_atoms);
    int n_types = header.n_types;

    // cell
    Cell cell(header.box, int_to_cell_type(header.cell_type));
    cfg.box_size = cell.matrix[0]; // approximate for orthorhombic

    // masses & charges
    std::vector<double> masses(n_types > 0 ? n_types : 1, 0.0);
    std::vector<double> charges(n_types > 0 ? n_types : 1, 0.0);
    if (n_types > 0) {
        is.read(reinterpret_cast<char*>(masses.data()), n_types * sizeof(double));
        is.read(reinterpret_cast<char*>(charges.data()), n_types * sizeof(double));
    }

    // ---- FF section ----
    is.seekg(static_cast<std::streamoff>(header.offset_ff));

    auto read_int_vec = [&](std::vector<int>& v, uint32_t count) {
        v.resize(count);
        is.read(reinterpret_cast<char*>(v.data()), count * sizeof(int32_t));
    };
    auto read_dbl_vec = [&](std::vector<double>& v, uint32_t count) {
        v.resize(count);
        is.read(reinterpret_cast<char*>(v.data()), count * sizeof(double));
    };

    // Bonds
    uint32_t n_bonds;
    is.read(reinterpret_cast<char*>(&n_bonds), sizeof(n_bonds));
    if (n_bonds > 0) {
        HarmonicBondParams bp;
        bp.n = n_bonds;
        read_int_vec(bp.i, n_bonds);
        read_int_vec(bp.j, n_bonds);
        read_dbl_vec(bp.k, n_bonds);
        read_dbl_vec(bp.r0, n_bonds);
        cfg.bonds.push_back(std::move(bp));
    }

    // Angles
    uint32_t n_angles;
    is.read(reinterpret_cast<char*>(&n_angles), sizeof(n_angles));
    if (n_angles > 0) {
        HarmonicAngleParams ap;
        ap.n = n_angles;
        read_int_vec(ap.i, n_angles);
        read_int_vec(ap.j, n_angles);
        read_int_vec(ap.k, n_angles);
        read_dbl_vec(ap.k_theta, n_angles);
        read_dbl_vec(ap.theta0, n_angles);
        cfg.angles.push_back(std::move(ap));
    }

    // Dihedrals
    uint32_t n_dihedrals;
    is.read(reinterpret_cast<char*>(&n_dihedrals), sizeof(n_dihedrals));
    if (n_dihedrals > 0) {
        PeriodicDihedralParams dp;
        dp.n = n_dihedrals;
        read_int_vec(dp.i, n_dihedrals);
        read_int_vec(dp.j, n_dihedrals);
        read_int_vec(dp.k, n_dihedrals);
        read_int_vec(dp.l, n_dihedrals);
        read_int_vec(dp.periodicity, n_dihedrals);
        read_dbl_vec(dp.k_phi, n_dihedrals);
        read_dbl_vec(dp.phi0, n_dihedrals);
        cfg.dihedrals.push_back(std::move(dp));
    }

    // LJ
    uint32_t n_lj_types;
    is.read(reinterpret_cast<char*>(&n_lj_types), sizeof(n_lj_types));
    if (n_lj_types > 0) {
        cfg.lj.n_types = n_lj_types;
        cfg.lj.sigma.resize(n_lj_types * n_lj_types);
        cfg.lj.epsilon.resize(n_lj_types * n_lj_types);
        is.read(reinterpret_cast<char*>(cfg.lj.sigma.data()), n_lj_types * n_lj_types * sizeof(double));
        is.read(reinterpret_cast<char*>(cfg.lj.epsilon.data()), n_lj_types * n_lj_types * sizeof(double));
        cfg.use_lj = true;
    }

    // ---- Positions section ----
    is.seekg(static_cast<std::streamoff>(header.offset_pos));

    auto read_float_vec = [&](std::vector<double>& v, uint64_t count) {
        std::vector<float> tmp(count);
        is.read(reinterpret_cast<char*>(tmp.data()), count * sizeof(float));
        v.assign(tmp.begin(), tmp.end());
    };

    uint64_t n64 = n;
    cfg.pos_x.resize(n); cfg.pos_y.resize(n); cfg.pos_z.resize(n);
    read_float_vec(cfg.pos_x, n64);
    read_float_vec(cfg.pos_y, n64);
    read_float_vec(cfg.pos_z, n64);

    if (header.pos_flags & POS_HAS_VEL) {
        cfg.vel_x.resize(n); cfg.vel_y.resize(n); cfg.vel_z.resize(n);
        read_float_vec(cfg.vel_x, n64);
        read_float_vec(cfg.vel_y, n64);
        read_float_vec(cfg.vel_z, n64);
    }

    // atom type indices (always present)
    std::vector<int32_t> type_idx(n);
    is.read(reinterpret_cast<char*>(type_idx.data()), n * sizeof(int32_t));

    cfg.masses.resize(n);
    cfg.charges.resize(n);
    for (size_t i = 0; i < n; ++i) {
        int t = type_idx[i];
        cfg.masses[i] = (n_types > 0 && t >= 0 && t < n_types) ? masses[t] : 1.0;
        cfg.charges[i] = (n_types > 0 && t >= 0 && t < n_types) ? charges[t] : 0.0;
    }

    return cfg;
}

// ---- Writer ----

void write_dmdin(const std::string& path, const SimulationConfig& cfg) {
    std::ofstream os(path, std::ios::binary);
    if (!os) throw std::runtime_error("Cannot write " + path);

    size_t n = cfg.pos_x.size();
    int n_types = static_cast<int>(cfg.lj.n_types);
    if (n_types == 0) n_types = 1;
    bool has_vel = !cfg.vel_x.empty();

    DMDinHeader header;
    header.magic = DMDIN_MAGIC;
    header.version = DMDIN_VERSION;
    header.n_atoms = n;
    header.cell_type = 0; // orthorhombic
    header.n_types = static_cast<int>(cfg.lj.n_types > 0 ? cfg.lj.n_types : 1);
    header.pos_flags = has_vel ? POS_HAS_VEL : 0;
    header.reserved = 0;
    header.offset_ff = sizeof(header) + n_types * 2 * sizeof(double);
    header.offset_pos = 0; // filled after FF section is written
    header.box[0] = cfg.box_size; header.box[1] = 0; header.box[2] = 0;
    header.box[3] = 0; header.box[4] = cfg.box_size; header.box[5] = 0;
    header.box[6] = 0; header.box[7] = 0; header.box[8] = cfg.box_size;

    os.write(reinterpret_cast<const char*>(&header), sizeof(header));

    // masses & charges: per-type from per-atom data
    std::vector<double> masses(n_types, 0.0);
    std::vector<double> charges(n_types, 0.0);
    std::vector<int> per_type_count(n_types, 0);
    if (!cfg.masses.empty() && !cfg.atom_types.empty()) {
        for (size_t i = 0; i < n && i < cfg.masses.size(); ++i) {
            int t = cfg.atom_types[i];
            if (t >= 0 && t < n_types) {
                masses[t] += cfg.masses[i];
                charges[t] += cfg.charges[i];
                per_type_count[t]++;
            }
        }
        for (int t = 0; t < n_types; ++t) {
            if (per_type_count[t] > 0) {
                masses[t] /= per_type_count[t];
                charges[t] /= per_type_count[t];
            } else {
                masses[t] = cfg.masses.empty() ? 1.0 : cfg.masses[0];
            }
        }
    } else {
        masses.assign(n_types, cfg.masses.empty() ? 1.0 : cfg.masses[0]);
    }
    os.write(reinterpret_cast<const char*>(masses.data()), n_types * sizeof(double));
    os.write(reinterpret_cast<const char*>(charges.data()), n_types * sizeof(double));

    // ---- FF section ----
    auto write_uint32 = [&](uint32_t v) { os.write(reinterpret_cast<const char*>(&v), sizeof(v)); };
    auto write_int_vec = [&](const std::vector<int>& v) {
        os.write(reinterpret_cast<const char*>(v.data()), v.size() * sizeof(int32_t));
    };
    auto write_dbl_vec = [&](const std::vector<double>& v) {
        os.write(reinterpret_cast<const char*>(v.data()), v.size() * sizeof(double));
    };

    // flatten all bond params
    {
        uint32_t total = 0;
        for (auto& b : cfg.bonds) total += static_cast<uint32_t>(b.n);
        write_uint32(total);
        for (auto& b : cfg.bonds) {
            write_int_vec(b.i); write_int_vec(b.j);
            write_dbl_vec(b.k); write_dbl_vec(b.r0);
        }
    }

    {
        uint32_t total = 0;
        for (auto& a : cfg.angles) total += static_cast<uint32_t>(a.n);
        write_uint32(total);
        for (auto& a : cfg.angles) {
            write_int_vec(a.i); write_int_vec(a.j); write_int_vec(a.k);
            write_dbl_vec(a.k_theta); write_dbl_vec(a.theta0);
        }
    }

    {
        uint32_t total = 0;
        for (auto& d : cfg.dihedrals) total += static_cast<uint32_t>(d.n);
        write_uint32(total);
        for (auto& d : cfg.dihedrals) {
            write_int_vec(d.i); write_int_vec(d.j); write_int_vec(d.k); write_int_vec(d.l);
            write_int_vec(d.periodicity);
            write_dbl_vec(d.k_phi); write_dbl_vec(d.phi0);
        }
    }

    uint32_t lj_n = static_cast<uint32_t>(cfg.lj.n_types > 0 ? cfg.lj.n_types : 1);
    write_uint32(lj_n);
    if (cfg.lj.n_types > 0) {
        os.write(reinterpret_cast<const char*>(cfg.lj.sigma.data()), lj_n * lj_n * sizeof(double));
        os.write(reinterpret_cast<const char*>(cfg.lj.epsilon.data()), lj_n * lj_n * sizeof(double));
    } else {
        double def_sigma = 0.3405, def_eps = 0.997;
        os.write(reinterpret_cast<const char*>(&def_sigma), sizeof(double));
        os.write(reinterpret_cast<const char*>(&def_eps), sizeof(double));
    }

    // ---- Positions section ----
    header.offset_pos = static_cast<uint64_t>(os.tellp());

    auto write_float = [&](const std::vector<double>& v) {
        std::vector<float> tmp(v.begin(), v.end());
        os.write(reinterpret_cast<const char*>(tmp.data()), tmp.size() * sizeof(float));
    };
    write_float(cfg.pos_x);
    write_float(cfg.pos_y);
    write_float(cfg.pos_z);
    if (has_vel) {
        write_float(cfg.vel_x);
        write_float(cfg.vel_y);
        write_float(cfg.vel_z);
    }
    // atom type indices
    std::vector<int32_t> types(n, 0);
    if (!cfg.atom_types.empty()) {
        for (size_t i = 0; i < n; ++i) types[i] = static_cast<int32_t>(cfg.atom_types[i]);
    }
    os.write(reinterpret_cast<const char*>(types.data()), n * sizeof(int32_t));

    // patch offset_pos in header
    os.seekp(offsetof(DMDinHeader, offset_pos));
    os.write(reinterpret_cast<const char*>(&header.offset_pos), sizeof(header.offset_pos));
}
