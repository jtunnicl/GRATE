cmake_minimum_required(VERSION 3.1.0)
project(Grate LANGUAGES CXX)

# options to build CLI and GUI versions
option(BUILD_GUI "Build Qt5 GUI version" ON)
option(BUILD_CLI "Build command line version" OFF)

# enable testing (disable with -DBUILD_TESTING=OFF)
include(CTest)

# default to RELEASE build if not specified
if (NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()

# option to gather profiling data
option(ENABLE_PROFILING "Enable gprof profiling" OFF)
if (ENABLE_PROFILING)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")

    # disable some optimisations as recommended here: https://github.com/jrfonseca/gprof2dot
    if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-omit-frame-pointer -fno-inline-functions -fno-inline-functions-called-once -fno-optimize-sibling-calls")
    endif()
endif()

# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# Instruct CMake to run moc automatically when needed
set(CMAKE_AUTOMOC ON)

# Create code from a list of Qt designer ui files
set(CMAKE_AUTOUIC ON)

# common sources
set(CPP_SOURCES
    tinyxml2/tinyxml2.cpp
    tinyxml2_wrapper.cpp
    gratetime.cpp
    sed.cpp
    riverprofile.cpp
    hydro.cpp
    model.cpp
)

# GUI specific sources
set(GUI_SOURCES
    qcustomplot.cpp
    mainwindow.cpp
    main.cpp
)
set(UI_SOURCES
    RwaveWin.ui
)

# CLI specific sources
set(CLI_SOURCES
    cli.cpp
)

# build library from common sources
add_library(grate_common ${CPP_SOURCES})
target_include_directories(grate_common PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# GUI build
if (BUILD_GUI)
    # Find the Qt libraries
    find_package(Qt5Core REQUIRED)
    find_package(Qt5Widgets REQUIRED)
    find_package(Qt5PrintSupport REQUIRED)
    set(QT5_LIBRARIES Qt5::Core Qt5::Widgets Qt5::PrintSupport)

    # Tell CMake to create the GrateRip executable
    add_executable(GrateRip ${UI_SOURCES} ${GUI_SOURCES})

    # Link against Qt 5 libraries
    target_link_libraries(GrateRip grate_common ${QT5_LIBRARIES})
endif()

# CLI build
if (BUILD_CLI)
    # create CLI executable and link
    add_executable(GrateRipCLI ${CLI_SOURCES})
    target_link_libraries(GrateRipCLI grate_common)
endif()

# copy input files to build dir for testing
configure_file(Input_Rip1_equil_1938.dat ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(hydro_series.dat ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(test_out.xml ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

# add tests
add_subdirectory(tests)
