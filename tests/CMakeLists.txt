# build compare program
add_executable(compare compare.cpp)

if (BUILD_CLI)
    # add the test of GrateRipCLI
    set(TEST_RUN_DIR ${CMAKE_CURRENT_BINARY_DIR}/GrateRipCLIFullTest)
    add_test(
        NAME GrateRipCLIFullTest
        COMMAND ${CMAKE_COMMAND}
            -DTEST_RUN_DIR=${TEST_RUN_DIR}
            -DTEST_SRC_DIR=${CMAKE_CURRENT_SOURCE_DIR}
            -DTEST_BINARY=$<TARGET_FILE:GrateRipCLI>
            -DCOMPARE_BINARY=$<TARGET_FILE:compare>
            -P ${CMAKE_CURRENT_SOURCE_DIR}/run_test.cmake
    )
endif (BUILD_CLI)
