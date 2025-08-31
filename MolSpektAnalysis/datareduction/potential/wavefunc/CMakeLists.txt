cmake_minimum_required(VERSION 3.1.0)

project(Wavefuncs)

set(CMAKE_AUTOMOC ON)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

add_library(Wavefuncs STATIC
	include/CoupledSineWaveFunc.h
	include/CouplingFuncs.h
    CoupledSineWaveFunc.cpp
	CouplingFuncs.cpp
)

target_include_directories(Wavefuncs PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(Wavefuncs
	FitObject
	MdiChild
	Molecule
	Potential
	Residualfit
	Qt5::Widgets
)
