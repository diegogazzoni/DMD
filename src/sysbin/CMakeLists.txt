add_library(dmd_sysbin
    dmdin.cpp
)
target_include_directories(dmd_sysbin PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(dmd_sysbin PUBLIC dmd_core dmd_force nlohmann_json::nlohmann_json)
