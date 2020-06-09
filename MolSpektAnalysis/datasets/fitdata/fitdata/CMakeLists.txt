cmake_minimum_required(VERSION 3.1.0)

project(FitData)

set(CMAKE_AUTOMOC ON)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

add_library(FitData STATIC
    fitdata.cpp
	fitdata.h
)

target_link_libraries(FitData Qt5::Widgets)
