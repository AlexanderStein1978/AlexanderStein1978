cmake_minimum_required(VERSION 3.1.0)

project(Datensatz)

set(CMAKE_AUTOMOC ON)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

add_library(Datensatz STATIC
	include/marker.h
	include/datensatz.h
    datensatz.cpp
	element.h
)

target_include_directories(Datensatz PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(Datensatz Qt5::Widgets)
