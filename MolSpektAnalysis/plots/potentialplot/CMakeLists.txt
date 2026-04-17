cmake_minimum_required(VERSION 3.5.0)

project(PotentialPlot)

set(CMAKE_AUTOMOC ON)

find_package(Qt6 COMPONENTS Widgets REQUIRED)

add_library(PotentialPlot STATIC
	include/potentialplot.h
    potentialplot.cpp
)

target_include_directories(PotentialPlot PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(PotentialPlot
	DiagWindow
	MolSpektAnalysisMainWindow
	Potential
	Residualfit
	TableWindow
<<<<<<< HEAD
	Qt6::Widgets
=======
	Qt5::Widgets
>>>>>>> f91263c093dfe7215b3249af2e1113f12a7a6877
)
