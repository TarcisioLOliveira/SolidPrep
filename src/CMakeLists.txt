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
)

add_subdirectory(pathfinding)
add_subdirectory(sizing)
add_subdirectory(element)
add_subdirectory(material)
add_subdirectory(meshing)
add_subdirectory(finite_element)
add_subdirectory(topology_optimization)
add_subdirectory(utils)
add_subdirectory(optimization)
add_subdirectory(density_filter)
add_subdirectory(projection)
add_subdirectory(global_stiffness_matrix)
add_subdirectory(function)
add_subdirectory(optimizer)

add_executable(${PROJECT_NAME} ${SOURCES})

target_link_libraries(${PROJECT_NAME} PUBLIC sizing pathfinding topology_optimization finite_element meshing element material optimization density_filter projection global_stiffness_matrix function optimizer utils ${GMSH_LIB} ${OpenCASCADE_LIBRARIES} ${LAPACKE_LIBRARIES} cblas ${BLAS_LIBRARIES} ${LAPACKE_LIBRARIES} ${MPI_LIBRARIES} MUMPS::MUMPS Eigen3::Eigen Threads::Threads)

if(OpenMP_CXX_FOUND)
    target_link_libraries(${PROJECT_NAME} PUBLIC OpenMP::OpenMP_CXX OpenMP::OpenMP_C OpenMP::OpenMP_Fortran)
endif()

install(TARGETS ${PROJECT_NAME} RUNTIME DESTINATION .
                                LIBRARY DESTINATION .
                                ARCHIVE DESTINATION .)
