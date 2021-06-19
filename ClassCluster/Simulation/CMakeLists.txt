cmake_minimum_required(VERSION 3.1.0)

project(Simulation)

set(CMAKE_AUTOMOC ON)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

add_library(Simulation STATIC
	include/vector.h
    include/window.h
    include/Calculation.h
    include/particle.h
    Calculation.cpp
    deltaesortfunctor.h
	vector.cpp
    vsortfunctor.h
    window.cpp
)

target_include_directories(Simulation PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(Simulation
    ControlWindow
    Potential
    Templates
    Tools
    Qt5::Widgets
)
