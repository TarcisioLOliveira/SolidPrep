cmake_minimum_required(VERSION 3.9.0)
set(PROJECT_NAME SolidPrep)

project(${PROJECT_NAME})

set(VERSION_MAJOR "0")
set(VERSION_MINOR "1")
set(VERSION_PATCH "1")
set(VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH})

enable_language(CXX)
enable_language(C)
enable_language(Fortran)
set(CMAKE_CXX_STANDARD 17)

set(NO_OPENMP OFF CACHE BOOL "Compile without OpenMP support")
set(THREADS_PREFER_PTHREAD_FLAG ON)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

find_package(OpenCASCADE REQUIRED)
if(OpenCASCADE_FOUND)
    message("-- OpenCASCADE found. Version: ${OpenCASCADE_VERSION_MAJOR}.${OpenCASCADE_VERSION_MINOR}.${OpenCASCADE_MAINTENANCE_VERSION}.")
    message("-- OpenCASCADE include directory: ${OpenCASCADE_INCLUDE_DIR}")
    message("-- OpenCASCADE libraries: ${OpenCASCADE_LIBRARIES}")
endif()
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(LAPACKE REQUIRED)
find_package(CBLAS REQUIRED)
find_package(Threads REQUIRED)
find_library(GMSH_LIB gmsh REQUIRED)
find_path(GMSH_INC gmsh.h REQUIRED)
find_package(VTK)
find_package(MPI REQUIRED)
find_package(MUMPS REQUIRED)

if(NOT NO_OPENMP)
    find_package(OpenMP)
    if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -fopenmp-simd")
    endif()
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas")
endif()

if(WIN32 OR CYGWIN)
    list(APPEND GMSH_LIB winmm wsock32 ws2_32 psapi)
endif()

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -march=native -pipe -fomit-frame-pointer -DNDEBUG")
    set(CMAKE_CXX_FLAGS_DEBUG  "${CMAKE_CXX_FLAGS_DEBUG} -g -Og")# -fsanitize=address")
endif()

# Enable link time optimization if supported
include(CheckIPOSupported)
check_ipo_supported(RESULT result OUTPUT output)
if(result)
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
        list(APPEND CMAKE_CXX_COMPILE_OPTIONS_IPO "-flto=auto")
    endif()
else()
    message(WARNING "IPO/LTO is not supported: ${output}")
endif()

include_directories(
    ${PROJECT_SOURCE_DIR}/include
    SYSTEM
    ${PROJECT_SOURCE_DIR}/rapidjson/include
    ${OpenCASCADE_INCLUDE_DIR}
    ${LAPACK_INCLUDE_DIR}
    ${BLAS_INCLUDE_DIR}
    ${LAPACKE_INCLUDE_DIR}
    ${GMSH_INC}
    ${VTK_INCLUDE_DIR}
    ${MPI_INCLUDE_PATH}
    ${MUMPS_INCLUDE_DIR}
)
   
add_subdirectory(src)
