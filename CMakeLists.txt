cmake_minimum_required(VERSION 3.15)
project(stretch_engine VERSION 1.0.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Set module path for custom find modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# Find required packages
find_package(FFTW3 REQUIRED)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# Enable optimizations
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Set optimization flags
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -march=native")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
elseif(MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /O2")
endif()

# Library source files
set(LIB_SOURCES
    src/core/AudioStretchProcessor.cpp
    src/core/STNDecomposition.cpp
    src/core/NoiseMorphing.cpp
    src/core/PitchScaling.cpp
)

# Create library target
add_library(stretch_engine ${LIB_SOURCES})

# Set include directories for internal headers
target_include_directories(stretch_engine
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src/core
)

# Set include directories
target_include_directories(stretch_engine
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
    PRIVATE
        ${FFTW3_INCLUDE_DIRS}
        ${EIGEN3_INCLUDE_DIR}
)

# Link dependencies
target_link_libraries(stretch_engine
    PUBLIC
        Eigen3::Eigen
    PRIVATE
        ${FFTW3_LIBRARIES}
)

# Create CLI executable target
add_executable(stretch_engine_cli src/cli/stretch_engine.cpp)
target_link_libraries(stretch_engine_cli PRIVATE stretch_engine)

# Installation rules
install(TARGETS stretch_engine stretch_engine_cli
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)

install(FILES include/stretch_engine.h
    DESTINATION include
)

# Testing configuration
enable_testing()
add_subdirectory(tests)

# Examples directory
add_executable(generate_test_files examples/generate_test_file.cpp)

# Create examples directory during build
add_custom_command(
    TARGET generate_test_files
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/examples
)

# Run the test file generator after building
add_custom_command(
    TARGET generate_test_files
    POST_BUILD
    COMMAND ./generate_test_files
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

# Install example WAV files
install(DIRECTORY ${CMAKE_BINARY_DIR}/examples/
    DESTINATION share/stretch_engine/examples
    FILES_MATCHING PATTERN "*.wav"
)

# Create examples directory in source
add_custom_command(
    TARGET generate_test_files
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_SOURCE_DIR}/examples
)

# Package configuration
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/stretch_engine-config-version.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

configure_package_config_file(
    "cmake/stretch_engine-config.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/stretch_engine-config.cmake"
    INSTALL_DESTINATION lib/cmake/stretch_engine
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/stretch_engine-config.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/stretch_engine-config-version.cmake"
    DESTINATION lib/cmake/stretch_engine
)
