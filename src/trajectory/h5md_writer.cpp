#include "h5md_writer.h"
#include "core/system_data.h"
#include "core/cell.h"
#include <hdf5.h>
#include <vector>
#include <stdexcept>
#include <cstring>

static std::string h5_error_str() {
    // HDF5 error stack auto-prints to stderr; return generic msg
    return "HDF5 operation failed";
}

static hid_t create_group(hid_t loc, const char* name) {
    hid_t g = H5Gcreate2(loc, name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (g < 0) throw std::runtime_error(std::string("H5Gcreate failed: ") + name);
    return g;
}

// Create an extendible 1D dataset for time or step
static hid_t create_time_step_ds(hid_t group, const char* name, hid_t dtype) {
    hsize_t dims[1] = {0};
    hsize_t max_dims[1] = {H5S_UNLIMITED};
    hid_t space = H5Screate_simple(1, dims, max_dims);
    if (space < 0) throw std::runtime_error(std::string("H5Screate failed: ") + name);

    hid_t plist = H5Pcreate(H5P_DATASET_CREATE);
    hsize_t chunk[1] = {1024};
    H5Pset_chunk(plist, 1, chunk);

    hid_t ds = H5Dcreate2(group, name, dtype, space, H5P_DEFAULT, plist, H5P_DEFAULT);
    H5Pclose(plist);
    H5Sclose(space);
    if (ds < 0) throw std::runtime_error(std::string("H5Dcreate failed: ") + name);
    return ds;
}

// Create an extendible 2D or 3D dataset for position/velocity/box
static hid_t create_value_ds(hid_t group, const char* name, hid_t dtype,
                              int ndim, const hsize_t* max_dims_inner) {
    hsize_t dims[3] = {0};
    hsize_t max_dims[3];
    max_dims[0] = H5S_UNLIMITED;
    for (int i = 1; i < ndim; ++i) max_dims[i] = max_dims_inner[i-1];

    hid_t space = H5Screate_simple(ndim, dims, max_dims);
    if (space < 0) throw std::runtime_error(std::string("H5Screate failed: ") + name);

    hid_t plist = H5Pcreate(H5P_DATASET_CREATE);
    hsize_t chunk[3] = {1};
    for (int i = 1; i < ndim; ++i) chunk[i] = max_dims[i];
    H5Pset_chunk(plist, ndim, chunk);

    hid_t ds = H5Dcreate2(group, name, dtype, space, H5P_DEFAULT, plist, H5P_DEFAULT);
    H5Pclose(plist);
    H5Sclose(space);
    if (ds < 0) throw std::runtime_error(std::string("H5Dcreate failed: ") + name);
    return ds;
}

// Extend and write a slice of a 1D dataset
static void append_scalar(hid_t dataset, hid_t dtype,
                           hsize_t frame, const void* value) {
    hsize_t new_size[1] = {frame + 1};
    if (H5Dset_extent(dataset, new_size) < 0)
        throw std::runtime_error("H5Dset_extent failed");

    hid_t fspace = H5Dget_space(dataset);
    hsize_t start[1] = {frame};
    hsize_t count[1] = {1};
    H5Sselect_hyperslab(fspace, H5S_SELECT_SET, start, nullptr, count, nullptr);

    hid_t mspace = H5Screate_simple(1, count, nullptr);

    if (H5Dwrite(dataset, dtype, mspace, fspace, H5P_DEFAULT, value) < 0) {
        H5Sclose(mspace);
        H5Sclose(fspace);
        throw std::runtime_error("H5Dwrite failed (scalar)");
    }
    H5Sclose(mspace);
    H5Sclose(fspace);
}

// Extend and write a 2D or 3D slice (first dim = frame, rest = inner dims)
static void append_slice(hid_t dataset, hid_t dtype,
                          hsize_t frame, int ndim,
                          const hsize_t* inner_dims, const void* data) {
    hsize_t new_size[3] = {frame + 1};
    for (int i = 1; i < ndim; ++i) new_size[i] = inner_dims[i-1];

    if (H5Dset_extent(dataset, new_size) < 0)
        throw std::runtime_error("H5Dset_extent failed");

    hid_t fspace = H5Dget_space(dataset);

    hsize_t start[3] = {frame, 0, 0};
    hsize_t count[3] = {1};
    for (int i = 1; i < ndim; ++i) count[i] = inner_dims[i-1];

    H5Sselect_hyperslab(fspace, H5S_SELECT_SET, start, nullptr, count, nullptr);

    hid_t mspace = H5Screate_simple(ndim, count, nullptr);

    if (H5Dwrite(dataset, dtype, mspace, fspace, H5P_DEFAULT, data) < 0) {
        H5Sclose(mspace);
        H5Sclose(fspace);
        throw std::runtime_error("H5Dwrite failed (slice)");
    }
    H5Sclose(mspace);
    H5Sclose(fspace);
}

H5MDWriter::H5MDWriter(const std::string& path, size_t n_atoms, bool write_vel)
    : n_atoms_(n_atoms), write_vel_(write_vel)
{
    file_ = H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file_ < 0)
        throw std::runtime_error("H5Fcreate failed: " + path);

    // /h5md group + version
    hid_t h5md = create_group(file_, "h5md");
    {
        int version[2] = {1, 1};
        hsize_t vdims[1] = {2};
        hid_t vspace = H5Screate_simple(1, vdims, nullptr);
        hid_t vds = H5Dcreate2(h5md, "version", H5T_NATIVE_INT, vspace,
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if (vds < 0) throw std::runtime_error("H5Dcreate version failed");
        H5Dwrite(vds, H5T_NATIVE_INT, vspace, vspace, H5P_DEFAULT, version);
        H5Dclose(vds);
        H5Sclose(vspace);
    }
    H5Gclose(h5md);

    // /particles/atoms
    hid_t particles = create_group(file_, "particles");
    hid_t atoms = create_group(particles, "atoms");
    H5Gclose(particles);

    // /particles/atoms/position
    hid_t pos_grp = create_group(atoms, "position");
    {
        hsize_t inner[2] = {n_atoms_, 3};
        pos_value_ = create_value_ds(pos_grp, "value", H5T_NATIVE_DOUBLE, 3, inner);
        pos_time_ = create_time_step_ds(pos_grp, "time", H5T_NATIVE_DOUBLE);
        pos_step_ = create_time_step_ds(pos_grp, "step", H5T_NATIVE_INT);
    }
    H5Gclose(pos_grp);

    // /particles/atoms/velocity (optional)
    if (write_vel_) {
        hid_t vel_grp = create_group(atoms, "velocity");
        {
            hsize_t inner[2] = {n_atoms_, 3};
            vel_value_ = create_value_ds(vel_grp, "value", H5T_NATIVE_DOUBLE, 3, inner);
            vel_time_ = create_time_step_ds(vel_grp, "time", H5T_NATIVE_DOUBLE);
            vel_step_ = create_time_step_ds(vel_grp, "step", H5T_NATIVE_INT);
        }
        H5Gclose(vel_grp);
    }

    // /particles/atoms/box/edges
    hid_t box_grp = create_group(atoms, "box");
    hid_t edges_grp = create_group(box_grp, "edges");
    {
        hsize_t inner[2] = {3, 3};
        box_value_ = create_value_ds(edges_grp, "value", H5T_NATIVE_DOUBLE, 3, inner);
        box_time_ = create_time_step_ds(edges_grp, "time", H5T_NATIVE_DOUBLE);
        box_step_ = create_time_step_ds(edges_grp, "step", H5T_NATIVE_INT);
    }
    H5Gclose(edges_grp);
    H5Gclose(box_grp);
    H5Gclose(atoms);
}

H5MDWriter::~H5MDWriter() {
    if (pos_value_ >= 0) H5Dclose(pos_value_);
    if (pos_time_ >= 0) H5Dclose(pos_time_);
    if (pos_step_ >= 0) H5Dclose(pos_step_);
    if (write_vel_) {
        if (vel_value_ >= 0) H5Dclose(vel_value_);
        if (vel_time_ >= 0) H5Dclose(vel_time_);
        if (vel_step_ >= 0) H5Dclose(vel_step_);
    }
    if (box_value_ >= 0) H5Dclose(box_value_);
    if (box_time_ >= 0) H5Dclose(box_time_);
    if (box_step_ >= 0) H5Dclose(box_step_);
    if (file_ >= 0) H5Fclose(file_);
}

void H5MDWriter::write_frame(const SystemData& sys, const Cell& cell) {
    hsize_t frame = frame_count_;

    // --- position ---
    {
        // Interleave pos_x, pos_y, pos_z into a flat array [n_atoms * 3]
        std::vector<double> buf(n_atoms_ * 3);
        for (size_t i = 0; i < n_atoms_; ++i) {
            buf[i * 3 + 0] = sys.pos_x[i];
            buf[i * 3 + 1] = sys.pos_y[i];
            buf[i * 3 + 2] = sys.pos_z[i];
        }
        hsize_t inner[2] = {n_atoms_, 3};
        append_slice(pos_value_, H5T_NATIVE_DOUBLE, frame, 3, inner, buf.data());
    }

    double t = sys.time;
    int st = static_cast<int>(sys.step);
    append_scalar(pos_time_, H5T_NATIVE_DOUBLE, frame, &t);
    append_scalar(pos_step_, H5T_NATIVE_INT, frame, &st);

    // --- velocity ---
    if (write_vel_) {
        std::vector<double> buf(n_atoms_ * 3);
        for (size_t i = 0; i < n_atoms_; ++i) {
            buf[i * 3 + 0] = sys.vel_x[i];
            buf[i * 3 + 1] = sys.vel_y[i];
            buf[i * 3 + 2] = sys.vel_z[i];
        }
        hsize_t inner[2] = {n_atoms_, 3};
        append_slice(vel_value_, H5T_NATIVE_DOUBLE, frame, 3, inner, buf.data());
        append_scalar(vel_time_, H5T_NATIVE_DOUBLE, frame, &t);
        append_scalar(vel_step_, H5T_NATIVE_INT, frame, &st);
    }

    // --- box/edges ---
    {
        double box[9];
        std::memcpy(box, cell.matrix, sizeof(double) * 9);
        hsize_t inner[2] = {3, 3};
        append_slice(box_value_, H5T_NATIVE_DOUBLE, frame, 3, inner, box);
        append_scalar(box_time_, H5T_NATIVE_DOUBLE, frame, &t);
        append_scalar(box_step_, H5T_NATIVE_INT, frame, &st);
    }

    ++frame_count_;
}
