set(PROJECT_NAME SolidPrep)

set(SOURCES_FINITE_ELEMENT
    direct_solver.cpp
    gradient_descent.cpp
    PCG.cpp
    mumps_solver.cpp
)

add_library(finite_element ${SOURCES_FINITE_ELEMENT})
