set(PROJECT_NAME SolidPrep)

set(SOURCES_GLOBAL_STIFFNESS_MATRIX
    lapack_dense_symmetric_banded.cpp
    mumps_sparse_symmetric.cpp
    eigen_sparse_symmetric.cpp
    eigen_sparse_asymmetric.cpp
)

add_library(global_stiffness_matrix ${SOURCES_GLOBAL_STIFFNESS_MATRIX})
