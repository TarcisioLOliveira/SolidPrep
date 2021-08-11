set(PROJECT_NAME SolidPrep)

set(SOURCES
    main.cpp
    project_data.cpp
    force.cpp
    support.cpp
    element.cpp
    beam_graph.cpp
    cross_section.cpp
    ground_structure.cpp
    utils.cpp
    material.cpp
    meshing.cpp
    finite_element.cpp
    visualization.cpp
)

add_subdirectory(pathfinding)
add_subdirectory(sizing)
add_subdirectory(element)
add_subdirectory(material)
add_subdirectory(meshing)
add_subdirectory(finite_element)
add_subdirectory(topology_optimization)
add_subdirectory(utils)

add_executable(${PROJECT_NAME} ${SOURCES} ${SOURCES_PATHFINDING} ${SOURCES_SIZING} ${SOURCES_ELEMENT} ${SOURCES_MATERIAL} ${SOURCES_MESHING} ${SOURCES_FINITE_ELEMENT} ${SOURCES_TOPOLOGY_OPTIMIZATION} ${SOURCES_UTILS})

target_link_libraries(${PROJECT_NAME} sizing pathfinding topology_optimization finite_element meshing element material ${GMSH_LIB} ${OpenCASCADE_LIBRARIES} ${LAPACKE_LIBRARIES} ${LAPACK_LIBRARIES} cblas nlopt)

install(TARGETS ${PROJECT_NAME} RUNTIME DESTINATION .
                                LIBRARY DESTINATION .
                                ARCHIVE DESTINATION .)
