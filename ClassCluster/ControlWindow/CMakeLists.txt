cmake_minimum_required(VERSION 3.1.0)

project(ControlWindow)

set(CMAKE_AUTOMOC ON)

find_package(Qt6 COMPONENTS Widgets REQUIRED)

add_library(ControlWindow STATIC
    include/controlwindow.h
    controlwindow.cpp
    potcontrol.cpp
    potcontrol.h
)

target_include_directories(ControlWindow PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(ControlWindow
    ParticleWatch
    Potential
    Simulation
    Qt6::Widgets
)
