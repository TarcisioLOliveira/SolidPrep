set(PROJECT_NAME SolidPrep)

set(SOURCES_SIZING
    beam_sizing.cpp
    standard_sizing.cpp
)

add_library(sizing OBJECT ${SOURCES_SIZING})
