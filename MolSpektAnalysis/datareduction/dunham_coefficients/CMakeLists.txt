cmake_minimum_required(VERSION 3.1.0)

project(DunhamCoefficients)

set(CMAKE_AUTOMOC ON)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

add_library(DunhamCoefficients STATIC
    duntable.cpp
	duntable.h
)

target_link_libraries(DunhamCoefficients Qt5::Widgets)
