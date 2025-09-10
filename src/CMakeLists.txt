set(PROJECT_NAME SolidPrep)

set(SOURCES
    main.cpp
    project_data.cpp
    force.cpp
    support.cpp
    element.cpp
    beam_graph.cpp
    cross_section.cpp
    utils.cpp
    material.cpp
    meshing.cpp
    finite_element.cpp
    visualization.cpp
    geometry.cpp
    view_handler.cpp
    projection.cpp
    optimizer.cpp
    multimaterial.cpp
    global_stiffness_matrix.cpp
    spring.cpp
    curvature.cpp
    internal_loads.cpp
    solver_manager.cpp
    shape_handler.cpp
)

add_subdirectory(pathfinding)
add_subdirectory(sizing)
add_subdirectory(element)
add_subdirectory(material)
add_subdirectory(meshing)
add_subdirectory(finite_element)
add_subdirectory(utils)
add_subdirectory(optimization)
add_subdirectory(density_filter)
add_subdirectory(projection)
add_subdirectory(global_stiffness_matrix)
add_subdirectory(function)
add_subdirectory(optimizer)
add_subdirectory(general_global_matrix)
add_subdirectory(general_solver)
add_subdirectory(boundary_element)
add_subdirectory(contact_element)
add_subdirectory(field)
add_subdirectory(nonlinear_solver)
add_subdirectory(math)
add_subdirectory(simulation)
add_subdirectory(project_specification)

add_executable(${PROJECT_NAME} ${SOURCES})

target_link_libraries(${PROJECT_NAME} PUBLIC
    project_specification
    simulation
    math
    sizing
    pathfinding 
    finite_element
    nonlinear_solver 
    meshing
    element
    boundary_element
    contact_element
    material
    field
    density_filter
    projection
    global_stiffness_matrix
    density_based_function
    node_shape_based_function
    density_based_optimizer
    node_shape_based_optimizer
    optimization
    general_solver
    general_global_matrix
    utils
    spview
    ${GMSH_LIB}
    ${OpenCASCADE_LIBRARIES}
    ${LAPACKE_LIBRARIES}
    cblas
    ${BLAS_LIBRARIES}
    ${LAPACKE_LIBRARIES}
    ${PETSC_LINK_LIBRARIES}
    ${MPI_LIBRARIES}
    MUMPS::MUMPS
    Eigen3::Eigen
    Threads::Threads
    ${GSL_LIBRARIES}
    jsoncpp_lib
    ${Teem_LIBRARY_DIRS}/libteem.so # ${Teem_LIBRARIES} does not use the full path
)

if(OpenMP_CXX_FOUND)
    target_link_libraries(${PROJECT_NAME} PUBLIC OpenMP::OpenMP_CXX OpenMP::OpenMP_C OpenMP::OpenMP_Fortran)
endif()

#if(USE_CUDA)
#    target_link_libraries(${PROJECT_NAME} PUBLIC NVHPC::MATH NVHPC::HOSTUTILS NVHPC::CUDA)
#endif()

install(TARGETS ${PROJECT_NAME} RUNTIME DESTINATION .
                                LIBRARY DESTINATION .
                                ARCHIVE DESTINATION .)
