set(PROJECT_NAME SolidPrep)

set(SOURCES_TOPOLOGY_OPTIMIZATION
    minimal_volume.cpp
    minimal_compliance.cpp
)

add_library(topology_optimization ${SOURCES_TOPOLOGY_OPTIMIZATION})
