cmake_minimum_required(VERSION 3.1.0)

project(TermTable)

set(CMAKE_AUTOMOC ON)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

add_library(TermTable STATIC
	include/cttidialog.h
	include/perturbation.h
	include/termtable.h
	include/termdata.h
	include/measuredtermenergies.h
    cttidialog.cpp
	measuredtermenergies.cpp
	termdata.cpp
	termtable.cpp
)

target_include_directories(TermTable PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(TermTable
	Molecule
	MolSpektAnalysisMainWindow
	TableWindow
	Tools
	Qt5::Widgets
)
