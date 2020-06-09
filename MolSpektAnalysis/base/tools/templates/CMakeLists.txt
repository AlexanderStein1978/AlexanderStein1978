cmake_minimum_required(VERSION 3.1.0)

project(Templates)

set(CMAKE_AUTOMOC ON)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

add_library(Templates STATIC
    doublevector.cpp
	doublevector.h
	flexvector.h
	heapsort.h
	matrix.h
	velement.h
)

target_link_libraries(Templates Qt5::Widgets)
