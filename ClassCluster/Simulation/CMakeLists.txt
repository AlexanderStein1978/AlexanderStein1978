cmake_minimum_required(VERSION 3.1.0)

project(Simulation)

set(CMAKE_AUTOMOC ON)

find_package(Qt5 COMPONENTS Widgets Network REQUIRED)

add_library(Simulation STATIC
    include/vector.h
    include/window.h
    include/Calculation.h
    include/particle.h
    include/potstruct.h
    Calculation.cpp
    deltaesortfunctor.h
    vector.cpp
    vsortfunctor.h
    window.cpp
	yz_sortfunctor.h
	random_sortfunctor.h
)

target_include_directories(Simulation PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(Simulation
    Calc
    DiagWindow
    LineTable
    MdiChild
    MonteCarloSim
    Network
    ParticleWatch
    PotentialDefinition
    PotentialStructs
    PotentialTypes
    Potential
    Residualfit
    TableWindow
    Templates
    Tools
    Qt5::Network
    Qt5::Widgets
)
