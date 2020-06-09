cmake_minimum_required(VERSION 3.1.0)

project(FitObject)

set(CMAKE_AUTOMOC ON)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

add_library(FitObject STATIC
    fitobject.cpp
	fitobject.h
	gaussian.cpp
	gaussian.h
	lorentzian.cpp
	lorentzian.h
)

target_link_libraries(FitObject Qt5::Widgets)
