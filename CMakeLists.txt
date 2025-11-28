cmake_minimum_required(VERSION 3.16)
project(coincfinder LANGUAGES CXX)

option(COINCFINDER_BUILD_TESTS "Build unit/integration tests" OFF)
option(COINCFINDER_BUILD_PYTHON "Build pybind11 Python module" ON)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Core library shared by CLIs, tests, and Python bindings.
add_library(coincfinder_core STATIC
    src/Coincidences.cpp
    src/ReadCSV.cpp
    src/RollingSingles.cpp
)
target_include_directories(coincfinder_core PUBLIC include)
target_compile_features(coincfinder_core PUBLIC cxx_std_20)

find_package(OpenMP QUIET)
if(OpenMP_CXX_FOUND)
  target_link_libraries(coincfinder_core PUBLIC OpenMP::OpenMP_CXX)
endif()

# CLI tools
add_executable(CoincFinder src/CoincFinder.cpp)
target_link_libraries(CoincFinder PRIVATE coincfinder_core)

add_executable(CoincPairs src/CoincPairs.cpp)
target_link_libraries(CoincPairs PRIVATE coincfinder_core)

# Tests
if(COINCFINDER_BUILD_TESTS)
  include(CTest)
  add_executable(CoincFinderTests src/CoincFinderTests.cpp)
  target_link_libraries(CoincFinderTests PRIVATE coincfinder_core)
  add_test(NAME coincfinder_unit COMMAND CoincFinderTests)

  add_executable(TestRolling Testing/TestRolling.cpp)
  target_link_libraries(TestRolling PRIVATE coincfinder_core)
endif()

# Python bindings
if(COINCFINDER_BUILD_PYTHON)
  find_package(pybind11 CONFIG QUIET)
  if(NOT pybind11_FOUND)
    message(FATAL_ERROR "pybind11 not found. Install pybind11 and expose its CMake config package.")
  endif()
  pybind11_add_module(coincfinder_py src/bindings.cpp)
  target_link_libraries(coincfinder_py PRIVATE coincfinder_core)
endif()

# Copy helper Python scripts to the build output for quick use with the .pyd/.so
file(GLOB COINCFINDER_PY_SCRIPTS "${CMAKE_SOURCE_DIR}/*.py")
if(COINCFINDER_PY_SCRIPTS)
  if(WIN32)
    set(COINCFINDER_PY_DEST "${CMAKE_BINARY_DIR}/Release")
  else()
    set(COINCFINDER_PY_DEST "${CMAKE_BINARY_DIR}")
  endif()
  add_custom_target(copy_python_helpers ALL
    COMMAND ${CMAKE_COMMAND} -E make_directory "${COINCFINDER_PY_DEST}"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${COINCFINDER_PY_SCRIPTS} "${COINCFINDER_PY_DEST}"
    COMMENT "Copying Python helper scripts to ${COINCFINDER_PY_DEST}")
endif()
