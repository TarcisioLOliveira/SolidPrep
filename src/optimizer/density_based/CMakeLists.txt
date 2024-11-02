set(PROJECT_NAME SolidPrep)

set(SOURCES_DENSITY_BASED_OPTIMIZER
    mma.cpp
    newton.cpp
)

add_library(density_based_optimizer ${SOURCES_DENSITY_BASED_OPTIMIZER})
