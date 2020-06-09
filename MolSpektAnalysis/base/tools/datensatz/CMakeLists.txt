cmake_minimum_required(VERSION 3.1.0)

project(Datensatz)

set(CMAKE_AUTOMOC ON)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

add_library(Datensatz STATIC
    datensatz.cpp
	datensatz.h
	element.h
	marker.h
)

target_link_libraries(Datensatz Qt5::Widgets)
