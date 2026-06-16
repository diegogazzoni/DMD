add_library(dmd_trajectory
    h5md_writer.cpp
)

target_include_directories(dmd_trajectory PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(dmd_trajectory PUBLIC dmd_core ${HDF5_C_LIBRARIES})
target_include_directories(dmd_trajectory PUBLIC ${HDF5_C_INCLUDE_DIRS})
