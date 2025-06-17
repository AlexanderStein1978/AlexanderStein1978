cmake_minimum_required(VERSION 3.1.0)

project(MonteCarloSim)

set(CMAKE_AUTOMOC ON)

find_package(Qt6 COMPONENTS Widgets REQUIRED)

add_library(MonteCarloSim STATIC
	include/fitresult.h
	include/includepot.h
	include/mcfsettingsdialog.h
	include/montecarlosim.h
	include/MonteCarloSimDialog.h
	include/processview.h
    fieldwindow.cpp
	fieldwindow.h
	fitresult.cpp
	mcfsettingsdialog.cpp
	montecarlosim.cpp
	MonteCarloSimDialog.cpp
	processview.cpp
)

target_include_directories(MonteCarloSim PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(MonteCarloSim PRIVATE
	MolSpektAnalysisMainWindow
	TableWindow
	Qt6::Widgets
)
