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
)

add_subdirectory(pathfinding)
add_subdirectory(sizing)
add_subdirectory(element)
add_subdirectory(material)
add_subdirectory(meshing)

add_executable(${PROJECT_NAME} ${SOURCES} ${SOURCES_PATHFINDING} ${SOURCES_SIZING} ${SOURCES_ELEMENT})

target_link_libraries(${PROJECT_NAME} sizing pathfinding meshing element material ${GMSH_LIB} ${OpenCASCADE_LIBRARIES} ${LAPACKE_LIBRARIES} ${LAPACK_LIBRARIES} cblas)

install(TARGETS ${PROJECT_NAME} RUNTIME DESTINATION .
                                LIBRARY DESTINATION .
                                ARCHIVE DESTINATION .)
