cmake_minimum_required(VERSION 3.1.0)

project(SpectSimulation)

set(CMAKE_AUTOMOC ON)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

add_library(SpectSimulation STATIC
    lineprofilepoint.cpp
	lineprofilepoint.h
	pointwiselineprofile.cpp
	pointwiselineprofile.h
	simline.h
	simulationstateproperty.h
	spectsimulation.cpp
	spectsimulation.h
)

target_link_libraries(SpectSimulation Qt5::Widgets)
