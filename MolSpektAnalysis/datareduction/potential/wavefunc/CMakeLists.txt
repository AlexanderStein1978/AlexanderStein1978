cmake_minimum_required(VERSION 3.5.0)

project(Wavefuncs)

set(CMAKE_AUTOMOC ON)

find_package(Qt6 COMPONENTS Widgets REQUIRED)

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
<<<<<<< HEAD
	Qt6::Widgets
=======
	Qt5::Widgets
>>>>>>> f91263c093dfe7215b3249af2e1113f12a7a6877
)
