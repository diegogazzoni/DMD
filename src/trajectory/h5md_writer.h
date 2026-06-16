#pragma once

#include <hdf5.h>
#include <string>
#include <cstddef>

struct SystemData;
struct Cell;

class H5MDWriter {
public:
    H5MDWriter(const std::string& path, size_t n_atoms, bool write_vel);
    ~H5MDWriter();

    H5MDWriter(const H5MDWriter&) = delete;
    H5MDWriter& operator=(const H5MDWriter&) = delete;

    void write_frame(const SystemData& sys, const Cell& cell);
    size_t frame_count() const { return frame_count_; }

private:
    hid_t file_{-1};

    hid_t pos_value_{-1}, pos_time_{-1}, pos_step_{-1};
    hid_t vel_value_{-1}, vel_time_{-1}, vel_step_{-1};
    hid_t box_value_{-1}, box_time_{-1}, box_step_{-1};

    size_t n_atoms_;
    size_t frame_count_{0};
    bool write_vel_;
};
