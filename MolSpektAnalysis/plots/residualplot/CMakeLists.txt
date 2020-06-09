cmake_minimum_required(VERSION 3.1.0)

project(ResidualPlot)

set(CMAKE_AUTOMOC ON)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

add_library(ResidualPlot STATIC
	PerturbationInfoTab.cpp
	PerturbationInfoTab.h
	PerturbationInfoTable.cpp
	PerturbationInfoTable.h
    residualplot.cpp
	residualplot.h
)

target_link_libraries(ResidualPlot Qt5::Widgets)
