set(PROJECT_NAME SolidPrep)

set(SOURCES_OPTIMIZATION
    GCMMASolver.cpp
    MMASolver.cpp
)

add_library(optimization ${SOURCES_OPTIMIZATION})
