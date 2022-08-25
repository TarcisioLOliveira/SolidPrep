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

add_executable(${PROJECT_NAME} ${SOURCES} ${SOURCES_PATHFINDING} ${SOURCES_SIZING} ${SOURCES_ELEMENT} ${SOURCES_MATERIAL} ${SOURCES_MESHING} ${SOURCES_FINITE_ELEMENT} ${SOURCES_TOPOLOGY_OPTIMIZATION} ${SOURCES_UTILS} ${SOURCES_OPTIMIZATION})

target_link_libraries(${PROJECT_NAME} sizing pathfinding topology_optimization finite_element meshing element material utils optimization ${GMSH_LIB} ${OpenCASCADE_LIBRARIES} ${LAPACKE_LIBRARIES} ${LAPACK_LIBRARIES} cblas)

install(TARGETS ${PROJECT_NAME} RUNTIME DESTINATION .
                                LIBRARY DESTINATION .
                                ARCHIVE DESTINATION .)
