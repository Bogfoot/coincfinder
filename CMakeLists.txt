cmake_minimum_required(VERSION 3.16)
project(CoincFinder LANGUAGES CXX)

# Default to Release on single-config generators (Makefiles, Ninja).
if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type selected; defaulting to Release.")
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)
endif()

# Force multi-config generators (Visual Studio, Xcode) to only offer Release.
if(CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_CONFIGURATION_TYPES "Release" CACHE STRING "Configs" FORCE)
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

add_library(coincfinder STATIC ${COINCFINDER_SOURCES})
target_include_directories(coincfinder PUBLIC ${PROJECT_SOURCE_DIR}/include)

# Optional OpenMP acceleration; works on both MSVC and GCC/Clang.
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(coincfinder PUBLIC OpenMP::OpenMP_CXX)
    target_compile_definitions(coincfinder PUBLIC COINCFINDER_WITH_OPENMP=1)
endif()

add_executable(CoincFinder src/CoincFinder.cpp)
target_link_libraries(CoincFinder PRIVATE coincfinder)

# Copy helper Python scripts next to the built artifacts.
file(GLOB PY_SCRIPTS "${PROJECT_SOURCE_DIR}/*.py")
if(PY_SCRIPTS)
    add_custom_command(TARGET CoincFinder POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PY_SCRIPTS} ${CMAKE_BINARY_DIR}
        COMMENT "Copying Python helper scripts to build directory")
endif()
