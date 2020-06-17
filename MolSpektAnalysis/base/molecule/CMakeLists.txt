cmake_minimum_required(VERSION 3.1.0)

project(Molecule)

set(CMAKE_AUTOMOC ON)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

add_library(Molecule STATIC
	include/atom.h
	include/Combobox.h
	include/elstate.h
	include/isotab.h
	include/molecule.h
	include/Table.h
	include/Transition.h
	Combobox.cpp
    atom.cpp
	elstate.cpp
	isotab.cpp
	molecule.cpp
	Table.cpp
	TableItem.h
	TableItem.cpp
	Transition.cpp
)

target_include_directories(Molecule PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(Molecule
	MdiChild
	Tools
	Qt5::Widgets
)
