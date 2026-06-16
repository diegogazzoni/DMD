add_library(dmd_barostat
    barostat.cpp
    berendsen_barostat.cpp
    andersen_barostat.cpp
)
target_include_directories(dmd_barostat PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(dmd_barostat PUBLIC dmd_core dmd_thermostat)
