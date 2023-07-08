set(PROJECT_NAME SolidPrep)

set(SOURCES_GLOBAL_STIFFNESS_MATRIX
    lapack_dense_symmetric_banded.cpp
    mumps_sparse_symmetric.cpp
    eigen_sparse_symmetric.cpp
    eigen_sparse_asymmetric.cpp
    petsc_sparse_symmetric.cpp
)

add_library(global_stiffness_matrix ${SOURCES_GLOBAL_STIFFNESS_MATRIX})

# if(USE_CUDA)
#     set(SOURCES_GLOBAL_STIFFNESS_MATRIX_CUDA
#         cusolver_sparse_symmetric.cpp
#     )
#     add_library(global_stiffness_matrix_cuda ${SOURCES_GLOBAL_STIFFNESS_MATRIX_CUDA})
#     
#     target_include_directories(
#         global_stiffness_matrix_cuda
#         PRIVATE
#         ${NVHPC_CUDA_INCLUDE_DIR}
#         ${NVHPC_ROOT_DIR}/math_libs/${NVHPC_CUDA_VERSION}/include
#         ${NVHPC_ROOT_DIR}/compilers/include/
#     )
#     target_link_libraries(global_stiffness_matrix_cuda PRIVATE utils NVHPC::MATH NVHPC::HOSTUTILS NVHPC::CUDA NVHPC::CUDART)
# 
#     target_link_libraries(global_stiffness_matrix PUBLIC global_stiffness_matrix_cuda)
# endif()
