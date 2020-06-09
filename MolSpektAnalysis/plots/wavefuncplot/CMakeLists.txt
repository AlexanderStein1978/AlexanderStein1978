cmake_minimum_required(VERSION 3.1.0)

project(WaveFuncPlot)

set(CMAKE_AUTOMOC ON)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

add_library(WaveFuncPlot STATIC
	wavefuncplot.cpp
	wavefuncplot.h
)

target_link_libraries(WaveFuncPlot Qt5::Widgets)
