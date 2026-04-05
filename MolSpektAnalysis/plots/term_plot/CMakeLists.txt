cmake_minimum_required(VERSION 3.5.0)

project(TermPlot)

set(CMAKE_AUTOMOC ON)

find_package(Qt6 COMPONENTS Widgets REQUIRED)

add_library(TermPlot STATIC
	include/TermPlot.h
	TermPlot.cpp
	TPViewList.h
)

target_include_directories(TermPlot PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(TermPlot
	DiagWindow
	FitData
	Templates
	Tools
	Qt6::Widgets
)
