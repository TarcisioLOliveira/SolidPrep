set(PROJECT_NAME SolidPrep)

set(SOURCES_MESHING
    gmsh.cpp
    standard_beam_mesher.cpp
)

add_library(meshing ${SOURCES_MESHING})
