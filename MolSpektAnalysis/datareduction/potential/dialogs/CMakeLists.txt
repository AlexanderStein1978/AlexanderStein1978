cmake_minimum_required(VERSION 3.1.0)

project(PotentialDialogs)

set(CMAKE_AUTOMOC ON)

find_package(Qt6 COMPONENTS Widgets REQUIRED)

add_library(PotentialDialogs STATIC
	include/aipsidiag.h
	include/CMLRPotDialog.h
	include/CoefficientDialog.h
	include/CreateAnaPotSeriesFromMCSSplinePotSeriesDialog.h
	include/CSWFImportDialog.h
	include/FitAnaPotDialog.h
	include/ImprovePotSeriesDialog.h
	include/SFQSCalcDialog.h
    aipsidiag.cpp
	CMLRPotDialog.cpp
	CoefficientDialog.cpp
	CreateAnaPotSeriesFromMCSSplinePotSeriesDialog.cpp
	CSWFImportDialog.cpp
	FitAnaPotDialog.cpp
	ImprovePotSeriesDialog.cpp
	SFQSCalcDialog.cpp
)

target_include_directories(PotentialDialogs PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(PotentialDialogs PRIVATE
    MolSpektAnalysisMainWindow
	Tools
	Qt6::Widgets
)
