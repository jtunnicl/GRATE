# build compare program
add_executable(compare compare.cpp)

# test the compare executable
macro(add_compare_test TEST_NAME)
    set(REF_FILE ${CMAKE_CURRENT_SOURCE_DIR}/test_compare/${TEST_NAME}_ref.txt)
    set(NEW_FILE ${CMAKE_CURRENT_SOURCE_DIR}/test_compare/${TEST_NAME}_new.txt)
    add_test(
        NAME ${TEST_NAME}
        COMMAND compare ${REF_FILE} ${NEW_FILE}
    )
endmacro(add_compare_test)

add_compare_test(testcompare1)
add_compare_test(testcompare2)
add_compare_test(testcompare3)
add_compare_test(testcompare4)
add_compare_test(testcompare5)
add_compare_test(testcompare6)
set_tests_properties(
    testcompare2
    testcompare3
    testcompare4
    testcompare5
    testcompare6
    PROPERTIES WILL_FAIL TRUE
)

# test GrateTime object
add_executable(test_gratetime test_gratetime.cpp)
target_link_libraries(test_gratetime grate_common)
add_test(
    NAME GrateTime
    COMMAND test_gratetime
)

# test the CLI version
if (BUILD_CLI)
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
