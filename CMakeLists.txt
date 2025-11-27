cmake_minimum_required(VERSION 3.16)
project(CoincFinder LANGUAGES CXX)

# Default to Release for single-config generators
if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Flattened output directories
foreach(TYPE RUNTIME LIBRARY ARCHIVE)
    set(CMAKE_${TYPE}_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
    foreach(CONFIG Release RelWithDebInfo MinSizeRel Debug)
        set(CMAKE_${TYPE}_OUTPUT_DIRECTORY_${CONFIG} ${CMAKE_BINARY_DIR}/${CONFIG})
    endforeach()
endforeach()

set(COINCFINDER_SOURCES
    src/Coincidences.cpp
    src/ReadCSV.cpp
    src/RollingSingles.cpp
)

add_library(coincfinder_lib STATIC ${COINCFINDER_SOURCES})
set_target_properties(coincfinder_lib PROPERTIES OUTPUT_NAME coincfinder)
target_include_directories(coincfinder_lib PUBLIC ${PROJECT_SOURCE_DIR}/include)

# Optional OpenMP (static on MSVC to avoid missing DLLs)
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(coincfinder_lib PUBLIC OpenMP::OpenMP_CXX)
    target_compile_definitions(coincfinder_lib PUBLIC COINCFINDER_WITH_OPENMP=1)
    if(MSVC)
        target_link_options(coincfinder_lib PUBLIC "/openmp:experimental" "/MT")
    endif()
endif()

add_executable(CoincFinder src/CoincFinder.cpp)
target_link_libraries(CoincFinder PRIVATE coincfinder_lib)

add_executable(CoincPairs src/CoincPairs.cpp)
target_link_libraries(CoincPairs PRIVATE coincfinder_lib)

# Pybind11 module (auto-fetch if missing)
option(BUILD_PYTHON_BINDINGS "Build coincfinder Python module" ON)
if(BUILD_PYTHON_BINDINGS)
    find_package(pybind11 CONFIG QUIET)
    if(NOT pybind11_FOUND)
        include(FetchContent)
        message(STATUS "pybind11 not found; fetching with FetchContent")
        FetchContent_Declare(
            pybind11
            GIT_REPOSITORY https://github.com/pybind/pybind11.git
            GIT_TAG        v2.12.0
        )
        FetchContent_MakeAvailable(pybind11)
    endif()

    pybind11_add_module(coincfinder_py MODULE src/bindings.cpp)
    target_link_libraries(coincfinder_py PRIVATE coincfinder_lib)
    set_target_properties(coincfinder_py PROPERTIES
        OUTPUT_NAME "coincfinder"
        SUFFIX $<$<BOOL:WIN32>:.pyd>
    )
    # Copy module to build/ for easy import
    add_custom_command(TARGET coincfinder_py POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            $<TARGET_FILE:coincfinder_py>
            ${CMAKE_BINARY_DIR}/$<TARGET_FILE_NAME:coincfinder_py>
        COMMENT "Copying Python module to build/ for import")
endif()

# Copy helper Python scripts next to build dir
file(GLOB PY_SCRIPTS "${PROJECT_SOURCE_DIR}/*.py")
if(PY_SCRIPTS)
    add_custom_target(copy_python ALL
        COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PY_SCRIPTS} ${CMAKE_BINARY_DIR}
        COMMENT "Copying Python helper scripts to build directory")
    add_dependencies(copy_python CoincFinder)
endif()
