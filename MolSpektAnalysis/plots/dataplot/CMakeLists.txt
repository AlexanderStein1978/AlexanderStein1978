cmake_minimum_required(VERSION 3.5.0)

project(DataPlot)

set(CMAKE_AUTOMOC ON)

find_package(Qt6 COMPONENTS Widgets REQUIRED)

add_library(DataPlot STATIC
	include/dataplot.h
    dataplot.cpp
	sortfuncs.cpp
	sortfuncs.h
)

target_include_directories(DataPlot PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(DataPlot
	DiagWindow
	MolSpektAnalysisMainWindow
<<<<<<< HEAD
	Qt6::Widgets
=======
	Qt5::Widgets
>>>>>>> f91263c093dfe7215b3249af2e1113f12a7a6877
)
