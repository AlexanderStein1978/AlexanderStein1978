cmake_minimum_required(VERSION 3.1.0)

project(WaveFuncPlot)

set(CMAKE_AUTOMOC ON)

find_package(Qt6 COMPONENTS Widgets REQUIRED)

add_library(WaveFuncPlot STATIC
	include/wavefuncplot.h
	wavefuncplot.cpp
)

target_include_directories(WaveFuncPlot PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(WaveFuncPlot
	DiagWindow
	MolSpektAnalysisMainWindow
	MonteCarloSim
	PotentialStructs
	Qt6::Widgets
)
