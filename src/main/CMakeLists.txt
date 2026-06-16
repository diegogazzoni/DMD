add_executable(dmd main.cpp)
target_link_libraries(dmd PRIVATE dmd_sim dmd_sysbin dmd_core)
target_include_directories(dmd PRIVATE ${CMAKE_SOURCE_DIR}/src)
