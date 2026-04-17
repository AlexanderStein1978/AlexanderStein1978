cmake_minimum_required(VERSION 3.5.0)

project(ResidualPlot)

set(CMAKE_AUTOMOC ON)

find_package(Qt6 COMPONENTS Widgets REQUIRED)

add_library(ResidualPlot STATIC
	include/residualplot.h
	PerturbationInfoTab.cpp
	PerturbationInfoTab.h
	PerturbationInfoTable.cpp
	PerturbationInfoTable.h
    residualplot.cpp
)

target_include_directories(ResidualPlot PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(ResidualPlot
	DiagWindow
	MolSpektAnalysisMainWindow
	Residualfit
<<<<<<< HEAD
	Qt6::Widgets
=======
	Qt5::Widgets
>>>>>>> f91263c093dfe7215b3249af2e1113f12a7a6877
)
