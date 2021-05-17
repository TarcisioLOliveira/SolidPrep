set(PROJECT_NAME SolidPrep)

set(SOURCES_MATERIAL
    linear_elastic_isotropic.cpp
    linear_elastic_orthotropic.cpp
)

add_library(material ${SOURCES_MATERIAL})
