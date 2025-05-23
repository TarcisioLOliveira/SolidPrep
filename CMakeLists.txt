cmake_minimum_required(VERSION 3.9.0)
set(PROJECT_NAME SolidPrep)

project(${PROJECT_NAME})

set(VERSION_MAJOR "0")
set(VERSION_MINOR "2")
set(VERSION_PATCH "0")
set(VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH})

enable_language(CXX)
enable_language(C)
enable_language(Fortran)
set(CMAKE_CXX_STANDARD 20)

set(BUILD_TESTS OFF CACHE BOOL "Compile tests (requires Catch2 v3)")
set(NO_OPENMP OFF CACHE BOOL "Compile without OpenMP support")
set(THREADS_PREFER_PTHREAD_FLAG ON)
set(USE_BLIS OFF CACHE BOOL "Use BLIS as BLAS library")
set(BLAS_PATH /opt/blis-native CACHE STRING "Path to BLAS files")
# set(USE_CUDA OFF CACHE BOOL "Compile with CUDA support")
# 
# if(USE_CUDA)
#     enable_language(CUDA)
# endif()

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# PkgConfig
find_package(PkgConfig REQUIRED)

find_package(jsoncpp REQUIRED)
get_target_property(JSON_INC_PATH jsoncpp_lib INTERFACE_INCLUDE_DIRECTORIES)
find_package(OpenCASCADE REQUIRED)
if(OpenCASCADE_FOUND)
    message("-- OpenCASCADE found. Version: ${OpenCASCADE_VERSION_MAJOR}.${OpenCASCADE_VERSION_MINOR}.${OpenCASCADE_MAINTENANCE_VERSION}.")
    message("-- OpenCASCADE include directory: ${OpenCASCADE_INCLUDE_DIR}")
    message("-- OpenCASCADE libraries: ${OpenCASCADE_LIBRARIES}")
endif()
if(USE_BLIS)
    set(BLA_VENDOR FLAME)
    add_compile_definitions(USE_BLIS)
endif()
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(LAPACKE REQUIRED)
find_package(CBLAS REQUIRED)
find_package(Threads REQUIRED)
find_library(GMSH_LIB gmsh REQUIRED PATHS /lib64/ /lib64/mpich/lib /lib64/openmpi/lib)
find_path(GMSH_INC gmsh.h REQUIRED PATHS /usr/include/ /usr/include/mpich-x86_64/)
find_package(VTK)
find_package(MPI REQUIRED)
find_package(MUMPS REQUIRED)
find_package(Eigen3 REQUIRED NO_MODULE)
find_package(GSL REQUIRED)

# PETSc
if (PKG_CONFIG_FOUND)
    pkg_check_modules(PETSC PETSc REQUIRED)
endif()

# if (PETSC_FOUND)
#     list(APPEND COMPILE_OPTIONS ${PETSC_CFLAGS})
# 
#     include_directories(${PETSC_INCLUDE_DIRS})
# 
#     set(LINK_FLAGS "${LINK_FLAGS} ${PETSC_LDFLAGS}")
# 
#     list(APPEND LIBRARIES ${PETSC_LINK_LIBRARIES})
# 
#     set(CMAKE_REQUIRED_FLAGS ${PETSC_CFLAGS})
# 
#     set(CMAKE_REQUIRED_INCLUDES "${PETSC_INCLUDE_DIRS}")
# endif()

if(NOT NO_OPENMP)
    find_package(OpenMP)
    if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -fopenmp-simd")
    endif()
    add_compile_definitions(MMA_WITH_OPENMP)
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas")
    add_compile_definitions(NO_OPENMP)
endif()

if(WIN32 OR CYGWIN)
    list(APPEND GMSH_LIB winmm wsock32 ws2_32 psapi)
endif()

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -m64")
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

add_subdirectory(modules/SPView)

include_directories(
    ${PROJECT_SOURCE_DIR}/include
    SYSTEM
    ${PROJECT_SOURCE_DIR}/modules/SPView/include
    ${OpenCASCADE_INCLUDE_DIR}
    ${LAPACK_INCLUDE_DIR}
    ${BLAS_INCLUDE_DIR}
    ${LAPACKE_INCLUDE_DIR}
    ${GMSH_INC}
    ${VTK_INCLUDE_DIR}
    ${MPI_INCLUDE_PATH}
    ${MUMPS_INCLUDE_DIRS}
    ${Eigen3_INCLUDE_DIR}
    ${PETSC_INCLUDE_DIRS}
    ${GSL_INCLUDE_DIR}
    ${JSON_INC_PATH}
)

if(USE_CUDA)
    add_compile_definitions(USE_CUDA)
endif()

add_subdirectory(src)

if(BUILD_TESTS)
    include(CTest)
    find_package(Catch2 2 REQUIRED)
    set(CMAKE_CTEST_ARGUMENTS "--output-on-failure")

    add_subdirectory(test)
endif()
