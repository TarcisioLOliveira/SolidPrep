set(PROJECT_NAME SolidPrep)

set(SOURCES_MATERIAL
    linear_elastic_isotropic.cpp
    linear_elastic_orthotropic.cpp
    mandible.cpp
    linear_elastic_orthotropic_field.cpp
)

add_library(material OBJECT ${SOURCES_MATERIAL})
