add_library(dmd_integrate
    integrator.cpp
    constraints.cpp
)

target_include_directories(dmd_integrate PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(dmd_integrate PUBLIC dmd_core)
