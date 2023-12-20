set(PROJECT_NAME SolidPrep)

set(SOURCES_BOUNDARY_ELEMENT
    BTRI3.cpp
    #BTRI6.cpp
)

add_library(boundary_element ${SOURCES_BOUNDARY_ELEMENT})
