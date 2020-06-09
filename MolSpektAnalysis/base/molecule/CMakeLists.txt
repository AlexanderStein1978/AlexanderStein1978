cmake_minimum_required(VERSION 3.1.0)

project(Molecule)

set(CMAKE_AUTOMOC ON)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

add_library(Molecule STATIC
    atom.cpp
	atom.h
	elstate.cpp
	elstate.h
	isotab.cpp
	isotab.h
	molecule.cpp
	molecule.h
)

target_link_libraries(Molecule Qt5::Widgets)
