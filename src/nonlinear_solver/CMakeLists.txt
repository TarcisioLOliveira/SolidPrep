set(PROJECT_NAME SolidPrep)

set(SOURCES_NONLINEAR_SOLVER
    steepest_descent.cpp
)

add_library(nonlinear_solver ${SOURCES_NONLINEAR_SOLVER})
