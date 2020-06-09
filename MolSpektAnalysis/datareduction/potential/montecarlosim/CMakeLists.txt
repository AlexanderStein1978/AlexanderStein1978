cmake_minimum_required(VERSION 3.1.0)

project(MonteCarloSim)

set(CMAKE_AUTOMOC ON)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

add_library(MonteCarloSim STATIC
    fieldwindow.cpp
	fieldwindow.h
	fitresult.cpp
	fitresult.h
	includepot.h
	mcfsettingsdialog.cpp
	mcfsettingsdialog.h
	montecarlosim.cpp
	montecarlosim.h
	MonteCarloSimDialog.cpp
	MonteCarloSimDialog.h
	processview.cpp
	processview.h
)

target_link_libraries(MonteCarloSim Qt5::Widgets)
