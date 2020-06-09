cmake_minimum_required(VERSION 3.1.0)

project(PotentialPlot)

set(CMAKE_AUTOMOC ON)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

add_library(PotentialPlot STATIC
    potentialplot.cpp
	potentialplot.h
)

target_link_libraries(PotentialPlot Qt5::Widgets)
