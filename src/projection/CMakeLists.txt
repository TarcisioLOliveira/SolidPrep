set(PROJECT_NAME SolidPrep)

set(SOURCES_PROJECTION
    none.cpp
    threshold.cpp
    heaviside.cpp
)

add_library(projection OBJECT ${SOURCES_PROJECTION})
