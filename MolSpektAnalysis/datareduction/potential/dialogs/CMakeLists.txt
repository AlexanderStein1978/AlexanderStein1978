cmake_minimum_required(VERSION 3.1.0)

project(PotentialDialogs)

set(CMAKE_AUTOMOC ON)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

add_library(PotentialDialogs STATIC
    aipsidiag.cpp
	aipsidiag.h
	CMLRPotDialog.cpp
	CMLRPotDialog.h
	CoefficientDialog.cpp
	CoefficientDialog.h
	CreateAnaPotSeriesFromMCSSplinePotSeriesDialog.cpp
	CreateAnaPotSeriesFromMCSSplinePotSeriesDialog.h
	CSWFImportDialog.cpp
	CSWFImportDialog.h
	FitAnaPotDialog.cpp
	FitAnaPotDialog.h
	ImprovePotSeriesDialog.cpp
	ImprovePotSeriesDialog.h
	SFQSCalcDialog.cpp
	SFQSCalcDialog.h
)

target_link_libraries(PotentialDialogs Qt5::Widgets)
