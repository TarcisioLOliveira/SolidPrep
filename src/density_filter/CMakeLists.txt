set(PROJECT_NAME SolidPrep)

set(SOURCES_DENSITY_FILTER
    convolution.cpp
    helmholtz.cpp
)

add_library(density_filter ${SOURCES_DENSITY_FILTER})
