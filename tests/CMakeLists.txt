enable_testing()

# Find Google Test
find_package(GTest REQUIRED)

# Basic test executable
add_executable(basic_tests
    basic_tests.cpp
)

# Link against your library and GTest
target_link_libraries(basic_tests
    PRIVATE
    stretch_engine  # Your library name
    GTest::gtest
    GTest::gtest_main
)

# Register tests
add_test(NAME basic_tests COMMAND basic_tests)

# Set working directory for tests
set_tests_properties(basic_tests
    PROPERTIES
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

# Add test data directory
add_custom_command(TARGET basic_tests POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_CURRENT_SOURCE_DIR}/test_data
        ${CMAKE_CURRENT_BINARY_DIR}/test_data
)

find_package(GTest REQUIRED)
