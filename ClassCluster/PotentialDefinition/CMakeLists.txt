cmake_minimum_required(VERSION 3.1.0)

project(PotentialDefinition)

set(CMAKE_AUTOMOC ON)

find_package(Qt5 COMPONENTS Widgets Network REQUIRED)

add_library(PotentialDefinition STATIC
    include/potentialdefinerinputdata.h
    include/potentialdefiner.h
    potentialdefinerinputdata.cpp
    potentialdefiner.cpp
)

target_include_directories(PotentialDefinition PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(PotentialDefinition
    Simulation
)
