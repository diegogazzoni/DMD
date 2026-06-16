add_library(dmd_force
    force_component.cpp
    harmonic_bond.cpp
    harmonic_angle.cpp
    periodic_dihedral.cpp
    lennard_jones.cpp
    coulomb_direct.cpp
    coulomb_pme.cpp
    force_engine.cpp
)

target_include_directories(dmd_force PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(dmd_force PUBLIC dmd_core)
