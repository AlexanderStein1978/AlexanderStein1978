cmake_minimum_required(VERSION 3.1.0)

project(SpectSimulation)

set(CMAKE_AUTOMOC ON)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

add_library(SpectSimulation STATIC
	include/simulationstateproperty.h
	include/spectsimulation.h
	include/pointwiselineprofile.h
	include/lineprofilepoint.h
    lineprofilepoint.cpp
	pointwiselineprofile.cpp
	simline.h
	spectsimulation.cpp
)

target_include_directories(SpectSimulation PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(SpectSimulation
	DiagWindow
	Potential
	Templates
	Qt5::Widgets
)
