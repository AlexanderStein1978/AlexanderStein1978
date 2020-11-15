cmake_minimum_required(VERSION 3.1.0)

project(ParticleWatch)

set(CMAKE_AUTOMOC ON)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

add_library(ParticleWatch STATIC
    include/particlewatchtable.h
    include/watchpoint.h
    particlecontribution.cpp
    particlecontribution.h
    particlewatchtable.cpp
    watchpoint.cpp
    watchstep.cpp
    watchstep.h
)

target_include_directories(ParticleWatch PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(ParticleWatch
    Simulation
    TableWindow
    Qt5::Widgets
)
