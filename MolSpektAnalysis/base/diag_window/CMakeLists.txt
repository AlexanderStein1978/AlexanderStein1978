cmake_minimum_required(VERSION 3.1.0)

project(DiagWindow)

set(CMAKE_AUTOMOC ON)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

add_library(DiagWindow STATIC
    DiagWindow.cpp
	DiagWindow.h
	Picture.cpp
	Picture.h
)

target_link_libraries(DiagWindow Qt5::Widgets)
