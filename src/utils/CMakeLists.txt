set(PROJECT_NAME SolidPrep)

set(SOURCES_UTILS
    csr.cpp
    sparse_matrix.cpp
    D_operations.cpp
    LIN3B.cpp
    basis_tensor.cpp
)

add_library(utils ${SOURCES_UTILS})
