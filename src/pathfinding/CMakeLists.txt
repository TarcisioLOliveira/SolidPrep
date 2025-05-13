set(PROJECT_NAME SolidPrep)

set(SOURCES_PATHFINDING
    meshless_astar.cpp
    visibility_graph.cpp
)

add_library(pathfinding OBJECT ${SOURCES_PATHFINDING})
