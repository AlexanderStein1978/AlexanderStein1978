cmake_minimum_required(VERSION 3.1.0)

project(PotentialControls)

set(CMAKE_AUTOMOC ON)

find_package(Qt6 COMPONENTS Widgets REQUIRED)

add_library(PotentialControls STATIC
	include/CreateAnaPotSeriesControl.h
	include/ImprovePotSeriesControl.h
	include/MonteCarloSimControl.h
	include/SFQSCalcControl.h
	include/ThreadControl.h
    CreateAnaPotSeriesControl.cpp
	ImprovePotSeriesControl.cpp
	MonteCarloSimControl.cpp
	SFQSCalcControl.cpp
	ThreadControl.cpp
)

target_include_directories(PotentialControls PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(PotentialControls
	MolSpektAnalysisMainWindow
	TableWindow
	Qt6::Widgets
)
