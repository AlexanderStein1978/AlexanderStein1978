cmake_minimum_required(VERSION 3.1.0)

project(IntensityHistogram)

set(CMAKE_AUTOMOC ON)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

add_library(IntensityHistogram STATIC
	include/intensityhistogram.h
    intensityhistogram.cpp
)

target_include_directories(IntensityHistogram PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(IntensityHistogram
	DiagWindow
	Spektrum
	Qt5::Widgets
)
