cmake_minimum_required(VERSION 3.1.0)

project(ResidualPlot)

set(CMAKE_AUTOMOC ON)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

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
	Qt5::Widgets
)
