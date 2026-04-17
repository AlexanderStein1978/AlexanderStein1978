cmake_minimum_required(VERSION 3.5.0)

project(Simulation)

set(CMAKE_AUTOMOC ON)

<<<<<<< HEAD
find_package(Qt6 COMPONENTS Widgets Network REQUIRED)
=======
find_package(Qt5 COMPONENTS Widgets Network REQUIRED)
>>>>>>> f91263c093dfe7215b3249af2e1113f12a7a6877

add_library(Simulation STATIC
    include/vector.h
    include/window.h
    include/Calculation.h
    include/particle.h
    include/potstruct.h
    Calculation.cpp
    vector.cpp
    vsortfunctor.h
    window.cpp
    potstruct.cpp
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
<<<<<<< HEAD
    Qt6::Network
    Qt6::Widgets
=======
    Qt5::Network
    Qt5::Widgets
>>>>>>> f91263c093dfe7215b3249af2e1113f12a7a6877
)
