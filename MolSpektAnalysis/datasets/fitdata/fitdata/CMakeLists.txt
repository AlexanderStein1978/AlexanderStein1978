cmake_minimum_required(VERSION 3.5.0)

project(FitData)

set(CMAKE_AUTOMOC ON)

find_package(Qt6 COMPONENTS Widgets REQUIRED)

add_library(FitData STATIC
	include/fitdatabasedata.h
	include/fitdata.h
	include/fitdatacore.h
    fitdata.cpp
    fitdatacore.cpp
)

target_include_directories(FitData PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(FitData
	FitDataDialogs
	FitDataSortFuncs
	FitDataStructs
	TableWindow
	Qt6::Widgets
)
