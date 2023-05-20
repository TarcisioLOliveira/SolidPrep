set(PROJECT_NAME SolidPrep)

set(SOURCES_PROJECTION
    none.cpp
    threshold.cpp
    heaviside.cpp
)

add_library(projection ${SOURCES_PROJECTION})
