cmake_minimum_required(VERSION 3.16)
project(CoincFinder LANGUAGES CXX)

# Default to Release on single-config generators (Makefiles, Ninja).
if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type selected; defaulting to Release.")
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)
endif()

# Prefer Release on multi-config generators; keep other configs available.
if(CMAKE_CONFIGURATION_TYPES)
    list(INSERT CMAKE_CONFIGURATION_TYPES 0 "Release")
    list(REMOVE_DUPLICATES CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_CONFIGURATION_TYPES "${CMAKE_CONFIGURATION_TYPES}" CACHE STRING "Configs" FORCE)
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Keep all outputs (libs + exe) in the build tree root for Windows & Linux.
foreach(TYPE RUNTIME LIBRARY ARCHIVE)
    set(CMAKE_${TYPE}_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
    foreach(CONFIG Release RelWithDebInfo MinSizeRel Debug)
        set(CMAKE_${TYPE}_OUTPUT_DIRECTORY_${CONFIG} ${CMAKE_BINARY_DIR})
    endforeach()
endforeach()

set(COINCFINDER_SOURCES
    src/Coincidences.cpp
    src/ReadCSV.cpp
)

# Use a name that is case-unique on Windows file systems.
add_library(coincfinder_lib STATIC ${COINCFINDER_SOURCES})
set_target_properties(coincfinder_lib PROPERTIES OUTPUT_NAME coincfinder)
target_include_directories(coincfinder_lib PUBLIC ${PROJECT_SOURCE_DIR}/include)

# Optional OpenMP acceleration; works on both MSVC and GCC/Clang.
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(coincfinder_lib PUBLIC OpenMP::OpenMP_CXX)
    target_compile_definitions(coincfinder_lib PUBLIC COINCFINDER_WITH_OPENMP=1)
endif()

add_executable(CoincFinder src/CoincFinder.cpp)
target_link_libraries(CoincFinder PRIVATE coincfinder_lib)

# Copy helper Python scripts next to the built artifacts after building the exe.
file(GLOB PY_SCRIPTS "${PROJECT_SOURCE_DIR}/*.py")
if(PY_SCRIPTS)
    add_custom_target(copy_python ALL
        COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PY_SCRIPTS} ${CMAKE_BINARY_DIR}
        COMMENT "Copying Python helper scripts to build directory")
    add_dependencies(copy_python CoincFinder)
endif()
