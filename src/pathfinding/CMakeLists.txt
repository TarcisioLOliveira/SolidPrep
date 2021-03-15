set(PROJECT_NAME SolidPrep)

set(SOURCES_PATHFINDING
    meshless_astar.cpp
)

add_library(pathfinding ${SOURCES_PATHFINDING})
