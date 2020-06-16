cmake_minimum_required(VERSION 3.1.0)

project(Potential)

set(CMAKE_AUTOMOC ON)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

add_library(Potential STATIC
	include/potential.h
	include/potworker.h
	include/PotFit.h
    potential.cpp
	PotFit.cpp
	PotFitLogWindow.cpp
	PotFitLogWindow.h
	potworker.cpp
)

target_include_directories(Potential PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(Potential
	MdiChild
	MonteCarloSim
	PotentialDialogs
	PotentialStructs
	PotentialTypes
	Qt5::Widgets
)
