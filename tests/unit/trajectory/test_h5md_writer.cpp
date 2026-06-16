#include "trajectory/h5md_writer.h"
#include "core/system_data.h"
#include "core/cell.h"
#include <gtest/gtest.h>
#include <hdf5.h>
#include <filesystem>

TEST(H5MDWriterTest, CreateAndWriteOneFrame) {
    std::string path = "test_one_frame.h5";
    {
        size_t n = 4;
        SystemData sys(n);
        sys.pos_x = {1.0, 2.0, 3.0, 4.0};
        sys.pos_y = {5.0, 6.0, 7.0, 8.0};
        sys.pos_z = {9.0, 10.0, 11.0, 12.0};
        sys.vel_x = {0.1, 0.2, 0.3, 0.4};
        sys.vel_y = {0.5, 0.6, 0.7, 0.8};
        sys.vel_z = {0.9, 1.0, 1.1, 1.2};
        sys.time = 0.0;
        sys.step = 0;

        double box[9] = {3.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 3.0};
        Cell cell(box, CellType::orthorhombic);

        H5MDWriter writer(path, n, /*write_vel=*/true);
        writer.write_frame(sys, cell);
        EXPECT_EQ(writer.frame_count(), 1);
    }

    // Verify the file exists and read back using HDF5 C API
    ASSERT_TRUE(std::filesystem::exists(path));

    hid_t file = H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    ASSERT_GE(file, 0);

    // Check /h5md/version
    hid_t ds = H5Dopen2(file, "/h5md/version", H5P_DEFAULT);
    ASSERT_GE(ds, 0);
    int version[2];
    H5Dread(ds, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, version);
    EXPECT_EQ(version[0], 1);
    EXPECT_EQ(version[1], 1);
    H5Dclose(ds);

    // Check position values
    ds = H5Dopen2(file, "/particles/atoms/position/value", H5P_DEFAULT);
    ASSERT_GE(ds, 0);
    hid_t space = H5Dget_space(ds);
    hsize_t dims[3];
    int ndim = H5Sget_simple_extent_ndims(space);
    H5Sget_simple_extent_dims(space, dims, nullptr);
    EXPECT_EQ(ndim, 3);
    EXPECT_EQ(dims[0], 1);  // 1 frame
    EXPECT_EQ(dims[1], 4);  // 4 atoms
    EXPECT_EQ(dims[2], 3);  // 3 coords
    double pos[4 * 3];
    H5Dread(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, pos);
    EXPECT_DOUBLE_EQ(pos[0], 1.0);
    EXPECT_DOUBLE_EQ(pos[1], 5.0);
    EXPECT_DOUBLE_EQ(pos[2], 9.0);
    EXPECT_DOUBLE_EQ(pos[9], 4.0);  // 4th atom x
    H5Sclose(space);
    H5Dclose(ds);

    // Check step
    ds = H5Dopen2(file, "/particles/atoms/position/step", H5P_DEFAULT);
    ASSERT_GE(ds, 0);
    int step_val;
    H5Dread(ds, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &step_val);
    EXPECT_EQ(step_val, 0);
    H5Dclose(ds);

    // Check velocity present
    ds = H5Dopen2(file, "/particles/atoms/velocity/value", H5P_DEFAULT);
    ASSERT_GE(ds, 0);
    H5Dclose(ds);

    // Check box
    ds = H5Dopen2(file, "/particles/atoms/box/edges/value", H5P_DEFAULT);
    ASSERT_GE(ds, 0);
    double box_val[9];
    H5Dread(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, box_val);
    EXPECT_DOUBLE_EQ(box_val[0], 3.0);
    EXPECT_DOUBLE_EQ(box_val[4], 3.0);
    EXPECT_DOUBLE_EQ(box_val[8], 3.0);
    H5Dclose(ds);

    H5Fclose(file);
    std::filesystem::remove(path);
}

TEST(H5MDWriterTest, MultipleFrames) {
    std::string path = "test_multi_frame.h5";
    {
        size_t n = 2;
        SystemData sys(n);
        Cell cell(3.0, 3.0, 3.0);

        H5MDWriter writer(path, n, /*write_vel=*/false);

        for (int step = 0; step < 5; ++step) {
            sys.pos_x[0] = 1.0 + step;
            sys.pos_y[0] = 2.0;
            sys.pos_z[0] = 3.0;
            sys.pos_x[1] = 4.0;
            sys.pos_y[1] = 5.0;
            sys.pos_z[1] = 6.0;
            sys.time = step * 0.002;
            sys.step = static_cast<size_t>(step);
            writer.write_frame(sys, cell);
        }
        EXPECT_EQ(writer.frame_count(), 5);
    }

    // Verify frame count
    hid_t file = H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    ASSERT_GE(file, 0);
    hid_t ds = H5Dopen2(file, "/particles/atoms/position/value", H5P_DEFAULT);
    hid_t space = H5Dget_space(ds);
    hsize_t dims[3];
    H5Sget_simple_extent_dims(space, dims, nullptr);
    EXPECT_EQ(dims[0], 5);
    H5Sclose(space);
    H5Dclose(ds);

    // Verify no velocity group
    hid_t g = H5Gopen2(file, "/particles/atoms/velocity", H5P_DEFAULT);
    EXPECT_LT(g, 0);  // should not exist
    // H5Gopen failed — don't close

    H5Fclose(file);
    std::filesystem::remove(path);
}

TEST(H5MDWriterTest, NoVelocityGroup) {
    std::string path = "test_novel.h5";
    {
        size_t n = 1;
        SystemData sys(n);
        Cell cell(3.0, 3.0, 3.0);
        H5MDWriter writer(path, n, /*write_vel=*/false);
        writer.write_frame(sys, cell);
    }
    // Verify velocity group is absent
    hid_t file = H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    hid_t g = H5Gopen2(file, "/particles/atoms/velocity", H5P_DEFAULT);
    EXPECT_LT(g, 0);
    H5Fclose(file);
    std::filesystem::remove(path);
}
