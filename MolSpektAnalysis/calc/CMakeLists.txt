cmake_minimum_required(VERSION 3.1.0)

project(Calc)

set(CMAKE_AUTOMOC ON)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

add_library(Calc STATIC
    fit.cpp
	fit.h
)

target_link_libraries(Calc Qt5::Widgets)
