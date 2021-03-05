cmake_minimum_required(VERSION 3.5.0)
set(PROJECT_NAME SolidPrep)

project(${PROJECT_NAME})

set(VERSION_MAJOR "0")
set(VERSION_MINOR "0")
set(VERSION_PATCH "0")
set(VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH})

set(CMAKE_CXX_STANDARD 11)

find_package(OpenCASCADE REQUIRED)
if(OpenCASCADE_FOUND)
    message("-- OpenCASCADE found. Version: ${OpenCASCADE_VERSION_MAJOR}.${OpenCASCADE_VERSION_MINOR}.${OpenCASCADE_MAINTENANCE_VERSION}.")
    message("-- OpenCASCADE include directory: ${OpenCASCADE_INCLUDE_DIR}")
    message("-- OpenCASCADE libraries: ${OpenCASCADE_LIBRARIES}")
endif()

include_directories(
    ${PROJECT_SOURCE_DIR}/include
    SYSTEM
    ${OpenCASCADE_INCLUDE_DIR}
)
   
add_subdirectory(src)
