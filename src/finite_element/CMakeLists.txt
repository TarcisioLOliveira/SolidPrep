set(PROJECT_NAME SolidPrep)

set(SOURCES_FINITE_ELEMENT
    direct_solver.cpp
    gradient_descent.cpp
    PCG.cpp
    mumps_solver.cpp
    eigen_pcg.cpp
    petsc_pcg.cpp
)

add_library(finite_element ${SOURCES_FINITE_ELEMENT})

if(USE_CUDA)
    set(SOURCES_FINITE_ELEMENT_CUDA
        cusolver.cpp
        cusolver_impl.cpp
    )
    add_library(finite_element_cuda ${SOURCES_FINITE_ELEMENT_CUDA})
    
    target_include_directories(
        finite_element_cuda
        PRIVATE
        ${NVHPC_CUDA_INCLUDE_DIR}
        ${NVHPC_ROOT_DIR}/math_libs/${NVHPC_CUDA_VERSION}/include
        ${NVHPC_ROOT_DIR}/compilers/include/
    )
    target_link_libraries(finite_element_cuda PRIVATE global_stiffness_matrix_cuda utils NVHPC::MATH NVHPC::HOSTUTILS NVHPC::CUDA NVHPC::CUDART)

    target_link_libraries(finite_element PUBLIC finite_element_cuda)
endif()

