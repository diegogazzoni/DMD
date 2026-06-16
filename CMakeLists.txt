cmake_minimum_required(VERSION 3.25)
project(dmd VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# ----- Options -----
option(DMD_BUILD_TESTS "Build unit tests" ON)

# ----- Dependencies -----
if(DMD_BUILD_TESTS)
    include(FetchContent)
    FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG v1.15.2
    )
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(googletest)
endif()

# ----- Subdirectories -----
add_subdirectory(src)

if(DMD_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()
