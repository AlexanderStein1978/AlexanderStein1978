cmake_minimum_required(VERSION 3.1.0)

project(Progressions)

set(CMAKE_AUTOMOC ON)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

add_library(Progressions STATIC
	aprog.h
	intprog.cpp
	intprog.h
	progression.h
    progressions.cpp
	progressions.h
	sortprog.h
	testprog.h
)

target_link_libraries(Progressions Qt5::Widgets)
