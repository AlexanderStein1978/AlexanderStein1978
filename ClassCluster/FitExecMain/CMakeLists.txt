cmake_minimum_required(VERSION 3.1.0)

project(FitExecMain)

add_definitions(-DDATA_DIRECTORY="${CMAKE_SOURCE_DIR}/ClassCluster/Data")

set(CMAKE_AUTOMOC ON)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

add_executable(FitExecMain
    fitexeccontrol.h
    fitexeccontrol.cpp
    main.cpp
)

target_link_libraries(FitExecMain
    ControlWindow
    FitObject
	Logging
    MolSpektAnalysisMainWindow
    ParticleWatch
    Simulation
    TestHelper
    Qt5::Widgets
)
