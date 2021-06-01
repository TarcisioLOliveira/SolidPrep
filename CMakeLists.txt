cmake_minimum_required(VERSION 3.5.0)
set(PROJECT_NAME SolidPrep)

project(${PROJECT_NAME})

set(VERSION_MAJOR "0")
set(VERSION_MINOR "0")
set(VERSION_PATCH "0")
set(VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH})

enable_language(CXX)
enable_language(Fortran)
set(CMAKE_CXX_STANDARD 11)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

find_package(OpenCASCADE REQUIRED)
if(OpenCASCADE_FOUND)
    message("-- OpenCASCADE found. Version: ${OpenCASCADE_VERSION_MAJOR}.${OpenCASCADE_VERSION_MINOR}.${OpenCASCADE_MAINTENANCE_VERSION}.")
    message("-- OpenCASCADE include directory: ${OpenCASCADE_INCLUDE_DIR}")
    message("-- OpenCASCADE libraries: ${OpenCASCADE_LIBRARIES}")
endif()
find_package(LAPACK REQUIRED)
find_package(LAPACKE REQUIRED)
find_package(CBLAS REQUIRED)
find_library(GMSH_LIB gmsh REQUIRED)
find_path(GMSH_INC gmsh.h REQUIRED)

if(WIN32 OR CYGWIN)
    list(APPEND GMSH_LIB winmm wsock32 ws2_32 psapi)
endif()


include_directories(
    ${PROJECT_SOURCE_DIR}/include
    SYSTEM
    ${OpenCASCADE_INCLUDE_DIR}
    ${LAPACK_INCLUDE_DIR}
    ${LAPACKE_INCLUDE_DIR}
    ${GMSH_INC}
)
   
add_subdirectory(src)
