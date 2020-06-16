cmake_minimum_required(VERSION 3.1.0)

project(DataPlot)

set(CMAKE_AUTOMOC ON)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

add_library(DataPlot STATIC
	include/dataplot.h
    dataplot.cpp
	sortfuncs.cpp
	sortfuncs.h
)

target_include_directories(DataPlot PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(DataPlot
	DiagWindow
	Qt5::Widgets
)
