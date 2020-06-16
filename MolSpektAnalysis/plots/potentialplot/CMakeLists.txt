cmake_minimum_required(VERSION 3.1.0)

project(PotentialPlot)

set(CMAKE_AUTOMOC ON)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

add_library(PotentialPlot STATIC
	include/potentialplot.h
    potentialplot.cpp
)

target_include_directories(PotentialPlot PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(PotentialPlot
	DiagWindow
	Residualfit
	Qt5::Widgets
)
