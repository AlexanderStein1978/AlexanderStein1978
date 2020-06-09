cmake_minimum_required(VERSION 3.1.0)

project(Wavefuncs)

set(CMAKE_AUTOMOC ON)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

add_library(Wavefuncs STATIC
    CoupledSineWaveFunc.cpp
	CoupledSineWaveFunc.h
	CouplingFuncs.cpp
	CouplingFuncs.h
)

target_link_libraries(Wavefuncs Qt5::Widgets)
