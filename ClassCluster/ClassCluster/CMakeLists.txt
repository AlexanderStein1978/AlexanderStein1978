cmake_minimum_required(VERSION 3.1.0)

project(ClassCluster)

set(CMAKE_AUTOMOC ON)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

add_executable(ClassCluster
    main.cpp
    Calculation.cpp
    Calculation.h
    controlwindow.cpp
    controlwindow.h
    deltaesortfunctor.h
    particle.h
	potcontrol.cpp
	potcontrol.h
    potstruct.h
    window.cpp
    window.h
)

target_link_libraries(ClassCluster
    Potential
    Templates
	Tools
    Qt5::Widgets
)
