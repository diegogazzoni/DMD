add_library(dmd_thermostat
    thermostat.cpp
    andersen.cpp
    berendsen.cpp
    nose_hoover.cpp
)
target_include_directories(dmd_thermostat PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(dmd_thermostat PUBLIC dmd_core)
