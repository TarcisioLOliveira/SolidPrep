set(PROJECT_NAME SolidPrep)

set(SOURCES_FIELD
    orthotropic_flow.cpp
    principal_stress.cpp
)

add_library(field OBJECT ${SOURCES_FIELD})

