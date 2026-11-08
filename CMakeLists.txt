cmake_minimum_required(VERSION 3.15)
project(CombAero CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Optimization and warning flags
if(MSVC)
    add_compile_options(/W4 /WX /permissive- /utf-8)
    add_compile_definitions(_CRT_SECURE_NO_WARNINGS)
else()
    add_compile_options(-Wall -Wextra -Wpedantic -Werror)
endif()

# -----------------------------------------------------------------------------
# Core Library: combaero_core
# -----------------------------------------------------------------------------
file(GLOB_RECURSE CORE_SOURCES src/*.cpp)
list(FILTER CORE_SOURCES EXCLUDE REGEX "src/main\\.cpp$")

add_library(combaero_core STATIC ${CORE_SOURCES})
target_include_directories(combaero_core PUBLIC include)
set_target_properties(combaero_core PROPERTIES POSITION_INDEPENDENT_CODE ON)

# -----------------------------------------------------------------------------
# Google Test Integration
# -----------------------------------------------------------------------------
# OFF for the Python build (pyproject.toml): it needs no tests, and the fetch
# below made every pip install depend on GitHub's archive host (#1059).
option(COMBAERO_ENABLE_TESTS "Build the C++ unit tests (fetches googletest)" ON)

if(COMBAERO_ENABLE_TESTS)
enable_testing()

include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/v1.15.2.zip
)
# For Windows: Do not override the internal heap flags of MSVC
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

# -----------------------------------------------------------------------------
# Core Units & Regression Tests
# -----------------------------------------------------------------------------

# Add all .cpp files in tests/ as individual test executables
file(GLOB TEST_SOURCES tests/*.cpp)

foreach(test_source ${TEST_SOURCES})
    get_filename_component(test_name ${test_source} NAME_WE)
    add_executable(${test_name} ${test_source})
    target_link_libraries(${test_name} PRIVATE combaero_core GTest::gtest_main)
    # Repo root on the include path so tests can pull in golden data headers
    # that live outside include/ (e.g. validation/*/data/*.h) without
    # duplicating reference numbers into the test file itself.
    target_include_directories(${test_name} PRIVATE ${CMAKE_SOURCE_DIR})
    add_test(NAME ${test_name} COMMAND ${test_name})
endforeach()
endif()

# -----------------------------------------------------------------------------
# Optional: Python bindings (enabled explicitly or via scikit-build-core)
# -----------------------------------------------------------------------------
option(COMBAERO_BUILD_PYTHON "Build Python extension" OFF)

if (COMBAERO_BUILD_PYTHON)
    set(PYBIND11_FINDPYTHON ON)
    find_package(pybind11 CONFIG REQUIRED)

    pybind11_add_module(_core python/combaero/_core.cpp)
    target_link_libraries(_core PRIVATE combaero_core)
    target_include_directories(_core PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/include
    )

    install(TARGETS _core DESTINATION combaero)
endif()

# -----------------------------------------------------------------------------
# Optional: Examples
# -----------------------------------------------------------------------------

option(COMBAERO_BUILD_EXAMPLES "Build standalone C++ examples" ON)

if(COMBAERO_BUILD_EXAMPLES)
    set(EXAMPLE_NAMES
        acoustics_example
        combustion_example
        combustion_state_example
        compressible_example
        cooling_example
        humidair_example
        incompressible_example
        orifice_example
        stagnation_example
        thermo_example
    )

    foreach(example ${EXAMPLE_NAMES})
        if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/examples/${example}.cpp")
            add_executable(${example} examples/${example}.cpp)
            target_link_libraries(${example} PRIVATE combaero_core)

            # Copy data for examples that need it (if available)
            if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/thermo_transport_data.json")
                add_custom_command(TARGET ${example} POST_BUILD
                    COMMAND ${CMAKE_COMMAND} -E copy_if_different
                    ${CMAKE_CURRENT_SOURCE_DIR}/thermo_transport_data.json
                    ${CMAKE_CURRENT_BINARY_DIR}/thermo_transport_data.json
                )
            endif()

            # Register each example as a smoke test. The examples are the only
            # C++ code exercised end-to-end the way a user writes it, and
            # building them is not enough: three of them segfaulted for some
            # time in mix() because they assigned the public State::X member
            # directly instead of calling set_X(), leaving Y empty. They
            # compiled and linked cleanly the whole time.
            add_test(NAME ${example}_smoke COMMAND ${example})
        endif()
    endforeach()
endif()
