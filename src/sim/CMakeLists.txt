add_library(dmd_sim
    simulation_engine.cpp
    checkpoint.cpp
)
target_include_directories(dmd_sim PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(dmd_sim PUBLIC dmd_force dmd_integrate)
