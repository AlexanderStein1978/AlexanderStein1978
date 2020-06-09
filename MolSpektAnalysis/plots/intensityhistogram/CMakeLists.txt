cmake_minimum_required(VERSION 3.1.0)

project(IntensityHistogram)

set(CMAKE_AUTOMOC ON)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

add_library(IntensityHistogram STATIC
    intensityhistogram.cpp
	intensityhistogram.h
)

target_link_libraries(IntensityHistogram Qt5::Widgets)
