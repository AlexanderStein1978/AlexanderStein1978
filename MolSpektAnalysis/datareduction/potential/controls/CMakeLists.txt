cmake_minimum_required(VERSION 3.1.0)

project(PotentialControls)

set(CMAKE_AUTOMOC ON)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

add_library(PotentialControls STATIC
    CreateAnaPotSeriesControl.cpp
	CreateAnaPotSeriesControl.h
	ImprovePotSeriesControl.cpp
	ImprovePotSeriesControl.h
	MonteCarloSimControl.cpp
	MonteCarloSimControl.h
	SFQSCalcControl.cpp
	SFQSCalcControl.h
	ThreadControl.cpp
	ThreadControl.h
)

target_link_libraries(PotentialControls Qt5::Widgets)
