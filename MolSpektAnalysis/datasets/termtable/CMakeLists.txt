cmake_minimum_required(VERSION 3.1.0)

project(TermTable)

set(CMAKE_AUTOMOC ON)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

add_library(TermTable STATIC
    cttidialog.cpp
	cttidialog.h
	measuredtermenergies.cpp
	measuredtermenergies.h
	perturbation.h
	termdata.cpp
	termdata.h
	termtable.cpp
	termtable.h
)

target_link_libraries(TermTable Qt5::Widgets)
