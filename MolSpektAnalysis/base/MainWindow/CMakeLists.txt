cmake_minimum_required(VERSION 3.1.0)

project(MolSpektAnalysisMainWindow)

set(CMAKE_AUTOMOC ON)

find_package(Qt6 COMPONENTS PrintSupport Widgets REQUIRED)

add_library(MolSpektAnalysisMainWindow STATIC
	include/MainWindow.h
	include/inputdialog.h
        inputdialog.cpp
	MainWindow.cpp
	showaction.h
)

target_include_directories(MolSpektAnalysisMainWindow PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(MolSpektAnalysisMainWindow
	AboutMolSpektAnalysis
	AddSpektrum
	DataPlot
	DunhamCoefficients
	IntensityHistogram
	PotentialControls
	PotentialDialogs
	PotentialPlot
	ResidualPlot
	TermTable
	TermView
	Tools
	WaveFuncPlot
	Qt6::PrintSupport
	Qt6::Widgets
)
