set(PROJECT_NAME SolidPrep)

set(SOURCES_SIMULATION
    marginal_bone_loss.cpp
)

add_library(simulation OBJECT ${SOURCES_SIMULATION})
