set(PROJECT_NAME SolidPrep)

set(SOURCES_ELEMENT
    beam_linear_2D.cpp
    GT9.cpp
    TRI3.cpp
    Q4.cpp
    Q4S.cpp
    TET4.cpp
)

add_library(element ${SOURCES_ELEMENT})
