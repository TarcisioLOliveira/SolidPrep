set(PROJECT_NAME SolidPrep)

set(SOURCES_FUNCTION
    compliance.cpp
    volume.cpp
    global_stress_pnorm_normalized.cpp
    global_stress_pnorm.cpp
    omni_machining.cpp
    am_support.cpp
    mass.cpp
    mass_first_material.cpp
    global_stress_heaviside.cpp
    mechanostat.cpp
)

add_library(function ${SOURCES_FUNCTION})
