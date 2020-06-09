cmake_minimum_required(VERSION 3.1.0)

project(TermView)

set(CMAKE_AUTOMOC ON)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

add_library(TermView STATIC
	termview.cpp
	termview.h
)

target_link_libraries(TermView Qt5::Widgets)
