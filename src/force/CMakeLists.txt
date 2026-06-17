find_library(FFTW3_OMP_LIB fftw3_omp)
if(FFTW3_OMP_LIB)
    add_library(fftw3_omp SHARED IMPORTED)
    set_target_properties(fftw3_omp PROPERTIES IMPORTED_LOCATION "${FFTW3_OMP_LIB}")
endif()

# OpenMP support via Homebrew libomp
# Apple Clang requires -Xclang -fopenmp + explicit include/link paths
find_library(OMP_LIB omp PATHS /opt/homebrew/opt/libomp/lib NO_DEFAULT_PATH)

add_library(dmd_force
    force_component.cpp
    harmonic_bond.cpp
    harmonic_angle.cpp
    periodic_dihedral.cpp
    lennard_jones.cpp
    coulomb_direct.cpp
    coulomb_pme.cpp
    coulomb_exclusion.cpp
    force_engine.cpp
)

target_include_directories(dmd_force PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(dmd_force PUBLIC dmd_core PkgConfig::FFTW)
if(OMP_LIB)
    target_compile_options(dmd_force PRIVATE -Xclang -fopenmp)
    target_include_directories(dmd_force PRIVATE /opt/homebrew/opt/libomp/include)
    target_link_libraries(dmd_force PRIVATE "${OMP_LIB}")
endif()
if(TARGET fftw3_omp)
    target_link_libraries(dmd_force PRIVATE fftw3_omp)
endif()
