cmake_minimum_required(VERSION 3.1.0)

project(FitData)

set(CMAKE_AUTOMOC ON)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

add_library(FitData STATIC
	include/fitdata.h
    fitdata.cpp
)

target_include_directories(FitData PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(FitData
	FitDataDialogs
	FitDataSortFuncs
	FitDataStructs
	TableWindow
	Qt5::Widgets
)
