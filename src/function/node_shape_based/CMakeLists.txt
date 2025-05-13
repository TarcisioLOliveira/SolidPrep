set(PROJECT_NAME SolidPrep)

set(SOURCES_NODE_SHAPE_BASED_FUNCTION
    compliance.cpp
    volume.cpp
    global_stress_heaviside.cpp
    mechanostat.cpp
)

add_library(node_shape_based_function OBJECT ${SOURCES_NODE_SHAPE_BASED_FUNCTION})
