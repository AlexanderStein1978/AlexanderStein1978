cmake_minimum_required(VERSION 3.1.0)

project(ClassCluster)

set(CMAKE_AUTOMOC ON)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

add_executable(ClassCluster
    Calculation.cpp
	Calculation.h
    controlwindow.cpp
	controlwindow.h
	deltaesortfunctor.h
	heapsort.h
	main.cpp
	particle.h
	picture.cpp
	picture.h
	window.cpp
	window.h
)

target_link_libraries(ClassCluster Qt5::Widgets)
