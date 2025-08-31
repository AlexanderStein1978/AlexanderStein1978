cmake_minimum_required(VERSION 3.1.0)

project(Templates)

set(CMAKE_AUTOMOC ON)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

add_library(Templates STATIC
	include/doublevector.h
	include/flexvector.h
	include/heapsort.h
	include/matrix.h
	include/velement.h
    doublevector.cpp
)

target_include_directories(Templates PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(Templates Qt5::Widgets)
