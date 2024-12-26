#=============================================================================
#
# CMake configuration file for the MathUtils project
#
#=============================================================================
cmake_minimum_required(VERSION 3.13)

project(MathUtils LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)

#=============================================================================
# Options
#=============================================================================
cmake_policy(SET CMP0077 NEW)

option(MATHUTILS_BUILD_TESTS "Build the MathUtils tests" ON)
option(MATHUTILS_BUILD_BOOST_TESTS "Build the MathUtils tests requiring Boost" ON)
option(MATHUTILS_SKIP_LONG_TESTS "Skip tests longer than 1 sec" OFF)
option(MATHUTILS_TEST_COVERAGE "Build for coverage" OFF)

#=============================================================================
# Compilation settings
#=============================================================================
if (${MATHUTILS_TEST_COVERAGE})
    add_compile_options("--coverage" "-O0")
    add_link_options("--coverage")
    add_custom_target(coverage gcovr . --root ../src --exclude-unreachable-branches --print-summary --html cov_report/ --html-details --xml coverage.xml --xml-pretty -j 4)
endif()

#=============================================================================
# Retrieving the current Git revision
#=============================================================================
set(PRE_CONFIGURE_FILE ${CMAKE_CURRENT_SOURCE_DIR}/src/MathUtils/version.h.in CACHE FILEPATH "" FORCE)
set(POST_CONFIGURE_FILE ${CMAKE_CURRENT_BINARY_DIR}/src/MathUtils/version.h CACHE FILEPATH "" FORCE)
set(GIT_WORKING_DIR ${CMAKE_CURRENT_SOURCE_DIR} CACHE FILEPATH "" FORCE)
set(GIT_STATE_FILE ${CMAKE_CURRENT_BINARY_DIR}/git-state-hash CACHE FILEPATH "" FORCE)
include(cmake/git_watcher.cmake)

if (NOT EXISTS ${POST_CONFIGURE_FILE})
    configure_file(${PRE_CONFIGURE_FILE} ${POST_CONFIGURE_FILE})
endif()

# =============================================================================
# MathUtils DEPENDENCIES
# =============================================================================
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/URL.conf.cmake)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

include(Add_Eigen3)
include(Add_Units)

#=============================================================================
# MathUtils library
#=============================================================================
add_subdirectory(src/MathUtils)

#=============================================================================
# Tests
#=============================================================================
if (${MATHUTILS_BUILD_TESTS})
    enable_testing()
    include(GoogleTest)
    include(Add_googletest)
    add_subdirectory(tests)
endif()
