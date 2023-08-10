set(PROJECT_NAME SolidPrep)

set(SOURCES_UTILS
    csr.cpp
    sparse_matrix.cpp
    D_operations.cpp
)

add_library(utils ${SOURCES_UTILS})
