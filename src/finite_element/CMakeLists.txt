set(PROJECT_NAME SolidPrep)

set(SOURCES_FINITE_ELEMENT
    direct_solver.cpp
    gradient_descent.cpp
)

add_library(finite_element ${SOURCES_FINITE_ELEMENT})
