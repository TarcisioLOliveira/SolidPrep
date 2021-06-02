set(PROJECT_NAME SolidPrep)

set(SOURCES_MESHING
    gmsh.cpp
)

add_library(meshing ${SOURCES_MESHING})
