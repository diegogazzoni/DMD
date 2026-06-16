add_library(dmd_core
    cell.cpp
    system_data.cpp
    error.cpp
)

target_include_directories(dmd_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/src)
