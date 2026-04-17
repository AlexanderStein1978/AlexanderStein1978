cmake_minimum_required(VERSION 3.5.0)

project(Potential)

set(CMAKE_AUTOMOC ON)

find_package(Qt6 COMPONENTS Widgets REQUIRED)

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

target_link_libraries(Potential PRIVATE
	Calc
	LineTable
	MdiChild
	MonteCarloSim
	PotentialDialogs
	PotentialStructs
	PotentialTypes
	TableWindow
	TermTable
	Tools
<<<<<<< HEAD
	Qt6::Widgets
=======
	Qt5::Widgets
>>>>>>> f91263c093dfe7215b3249af2e1113f12a7a6877
)
