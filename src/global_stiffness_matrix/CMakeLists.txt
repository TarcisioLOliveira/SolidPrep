set(PROJECT_NAME SolidPrep)

set(SOURCES_GLOBAL_STIFFNESS_MATRIX
    lapack_dense_symmetric_banded.cpp
)

add_library(global_stiffness_matrix ${SOURCES_GLOBAL_STIFFNESS_MATRIX})
