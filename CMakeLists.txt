cmake_minimum_required(VERSION 3.1.0)
project(Grate LANGUAGES CXX)

# option to gather profiling data
option(ENABLE_PROFILING "Enable gprof profiling" OFF)
if (ENABLE_PROFILING)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
endif()

# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# Instruct CMake to run moc automatically when needed
set(CMAKE_AUTOMOC ON)

# Create code from a list of Qt designer ui files
set(CMAKE_AUTOUIC ON)

# Find the Qt libraries
find_package(Qt5Widgets REQUIRED)
find_package(Qt5PrintSupport REQUIRED)
set(QT5_LIBRARIES Qt5::Widgets Qt5::PrintSupport)

# Populate CMake variables with the sources
set(CPP_SOURCES
    tinyxml2/tinyxml2.cpp
    sed.cpp
    riverprofile.cpp
    qcustomplot.cpp
    hydro.cpp
    mainwindow.cpp
    main.cpp
)
set(UI_SOURCES
    RwaveWin.ui
)

# Tell CMake to create the GrateRip executable
add_executable(GrateRip ${UI_SOURCES} ${CPP_SOURCES})

# Link against Qt 5 libraries
target_link_libraries(GrateRip ${QT5_LIBRARIES})

# copy input files to build dir for testing
configure_file(Input_Rip1_equil_1938.dat ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(hydro_series.dat ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(sed_series.dat ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(test_out.xml ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
