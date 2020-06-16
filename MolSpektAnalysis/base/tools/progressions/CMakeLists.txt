cmake_minimum_required(VERSION 3.1.0)

project(Progressions)

set(CMAKE_AUTOMOC ON)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

add_library(Progressions STATIC
	include/aprog.h
	include/intprog.h
	include/progression.h
	include/progressions.h
	include/sortprog.h
	include/testprog.h
	progressions.cpp
	intprog.cpp
)

target_include_directories(Progressions PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(Progressions
	Datensatz
	Molecule
	MolSpektAnalysisStructs
	Qt5::Widgets
)
