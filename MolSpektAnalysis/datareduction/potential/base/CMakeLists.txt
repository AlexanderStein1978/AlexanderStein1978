cmake_minimum_required(VERSION 3.1.0)

project(Potential)

set(CMAKE_AUTOMOC ON)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

add_library(Potential STATIC
    potential.cpp
	potential.h
	PotFit.cpp
	PotFit.h
	PotFitLogWindow.cpp
	PotFitLogWindow.h
	potworker.cpp
	potworker.h
)

target_link_libraries(Potential Qt5::Widgets)
