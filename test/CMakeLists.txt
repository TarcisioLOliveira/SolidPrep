set(TEST_EXEC_NAME ${PROJECT_NAME}_tests)


set(SOURCE_TESTS 
    main.cpp
    matrix.cpp
)

add_executable(${TEST_EXEC_NAME} ${SOURCE_TESTS})
target_link_libraries(${TEST_EXEC_NAME} PUBLIC
    Catch2::Catch2
    math
    ${LAPACKE_LIBRARIES}
    cblas
    ${BLAS_LIBRARIES}
    ${LAPACKE_LIBRARIES}
)

include(Catch)
catch_discover_tests(${TEST_EXEC_NAME})
