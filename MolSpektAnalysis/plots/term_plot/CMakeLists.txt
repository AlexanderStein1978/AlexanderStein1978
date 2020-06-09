cmake_minimum_required(VERSION 3.1.0)

project(TermPlot)

set(CMAKE_AUTOMOC ON)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

add_library(TermPlot STATIC
	TermPlot.cpp
	TermPlot.h
	TPViewList.h
)

target_link_libraries(TermPlot Qt5::Widgets)
