set(PROJECT_NAME SolidPrep)

set(SOURCES
    main.cpp
    project_data.cpp
    force.cpp
    support.cpp
)

add_subdirectory(pathfinding)

add_executable(${PROJECT_NAME} ${SOURCES} ${SOURCES_PATHFINDING})

target_link_libraries(${PROJECT_NAME} ${OpenCASCADE_LIBRARIES} pathfinding)

install(TARGETS ${PROJECT_NAME} RUNTIME DESTINATION .
                                LIBRARY DESTINATION .
                                ARCHIVE DESTINATION .)
