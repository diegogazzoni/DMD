add_executable(test_h5md_writer test_h5md_writer.cpp)
target_link_libraries(test_h5md_writer PRIVATE dmd_trajectory dmd_core ${HDF5_C_LIBRARIES} GTest::gtest GTest::gtest_main)
target_include_directories(test_h5md_writer PRIVATE ${HDF5_C_INCLUDE_DIRS})
add_test(NAME test_h5md_writer COMMAND test_h5md_writer)
