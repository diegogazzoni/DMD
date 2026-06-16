add_executable(test_cell test_cell.cpp)
target_link_libraries(test_cell PRIVATE dmd_core gtest_main)
add_test(NAME test_cell COMMAND test_cell)

add_executable(test_system_data test_system_data.cpp)
target_link_libraries(test_system_data PRIVATE dmd_core gtest_main)
add_test(NAME test_system_data COMMAND test_system_data)
