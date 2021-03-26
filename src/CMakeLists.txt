set(PROJECT_NAME SolidPrep)

set(SOURCES
    main.cpp
    project_data.cpp
    force.cpp
    support.cpp
    internal_force.cpp
)

add_subdirectory(pathfinding)
add_subdirectory(sizing)

add_executable(${PROJECT_NAME} ${SOURCES} ${SOURCES_PATHFINDING} ${SOURCES_SIZING})

SET(CMAKE_EXE_LINKER_FLAGS ${LAPACKE_LINKER_FLAGS})

target_link_libraries(${PROJECT_NAME} pathfinding sizing ${OpenCASCADE_LIBRARIES} ${LAPACKE_LIBRARIES} ${LAPACK_LIBRARIES})

install(TARGETS ${PROJECT_NAME} RUNTIME DESTINATION .
                                LIBRARY DESTINATION .
                                ARCHIVE DESTINATION .)
