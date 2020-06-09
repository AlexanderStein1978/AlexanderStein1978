cmake_minimum_required(VERSION 3.1.0)

project(MdiChild)

set(CMAKE_AUTOMOC ON)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

add_library(MdiChild STATIC
    mdichild.cpp
	mdichild.h
)

target_link_libraries(MdiChild Qt5::Widgets)
