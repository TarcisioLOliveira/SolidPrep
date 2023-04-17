set(PROJECT_NAME SolidPrep)

set(SOURCES_GENERAL_SOLVER
    petsc_general_pcg.cpp
)

add_library(general_solver ${SOURCES_GENERAL_SOLVER})
