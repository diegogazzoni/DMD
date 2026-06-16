add_executable(test_dmdin test_dmdin.cpp)
target_link_libraries(test_dmdin PRIVATE dmd_sysbin dmd_core gtest_main)
add_test(NAME test_dmdin COMMAND test_dmdin)
