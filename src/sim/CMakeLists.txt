add_library(dmd_sim
    simulation_engine.cpp
    checkpoint.cpp
    json_config.cpp
)
target_include_directories(dmd_sim PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(dmd_sim PUBLIC dmd_force dmd_integrate dmd_thermostat dmd_barostat dmd_trajectory nlohmann_json::nlohmann_json ${HDF5_C_LIBRARIES})
target_include_directories(dmd_sim PUBLIC ${HDF5_C_INCLUDE_DIRS})
