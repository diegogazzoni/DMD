# pybind11 bindings for DMD C++ core
# Builds _dmd_core.so which is imported by the Python 'dmd' package.
# Use: pip install -e .  or  cmake --build build && cp build/python/dmd/_dmd_core*.so python/dmd/

include(FetchContent)

FetchContent_Declare(
    pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11.git
    GIT_TAG v2.13.6
)
FetchContent_MakeAvailable(pybind11)

pybind11_add_module(_dmd_core MODULE
    dmd/core.cpp
)

target_include_directories(_dmd_core
    PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${pybind11_INCLUDE_DIRS}
)

target_link_libraries(_dmd_core
    PRIVATE
        dmd_sim
        dmd_core
        ${HDF5_C_LIBRARIES}
)

target_include_directories(_dmd_core
    PRIVATE
        ${HDF5_C_INCLUDE_DIRS}
)

# Copy the built module into the Python package directory
# so pip install -e . picks it up
set_target_properties(_dmd_core PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/python/dmd"
)
