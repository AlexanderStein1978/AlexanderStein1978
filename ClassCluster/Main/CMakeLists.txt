cmake_minimum_required(VERSION 3.1.0)

project(ClassCluster)

set(CMAKE_AUTOMOC ON)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

add_executable(ClassCluster
    main.cpp
)

target_link_libraries(ClassCluster
    ControlWindow
    FitObject
	Logging
    MolSpektAnalysisMainWindow
    ParticleWatch
    Simulation
    Qt5::Widgets
)
