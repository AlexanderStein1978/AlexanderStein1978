cmake_minimum_required(VERSION 3.1.0)

project(LineTable)

set(CMAKE_AUTOMOC ON)

find_package(Qt6 COMPONENTS Widgets REQUIRED)

add_library(LineTable STATIC
	include/fspdiag.h
	include/linetable.h
	include/linetablesortfunctions.h
	include/tableline.h
	include/termenergy.h
    fspdiag.cpp
	linetable.cpp
	questionbox.cpp
	questionbox.h
	linetablesortfunctions.cpp
	vsolistelement.h
)

target_include_directories(LineTable PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(LineTable PRIVATE
	Molecule
	MolSpektAnalysisMainWindow
	PotentialPlot
	Progressions
	Spektrum
	Qt6::Widgets
)
