set(PROJECT_NAME SolidPrep)

set(SOURCES_OPTIMIZER
    mma.cpp
    newton.cpp
)

add_library(optimizer ${SOURCES_OPTIMIZER})
