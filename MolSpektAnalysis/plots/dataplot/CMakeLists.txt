cmake_minimum_required(VERSION 3.1.0)

project(DataPlot)

set(CMAKE_AUTOMOC ON)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

add_library(DataPlot STATIC
    dataplot.cpp
	dataplot.h
	sortfuncs.cpp
	sortfuncs.h
)

target_link_libraries(DataPlot Qt5::Widgets)
