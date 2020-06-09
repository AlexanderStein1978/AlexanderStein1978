cmake_minimum_required(VERSION 3.1.0)

project(LineTable)

set(CMAKE_AUTOMOC ON)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

add_library(LineTable STATIC
    fspdiag.cpp
	fspdiag.h
	linetable.cpp
	linetable.h
	questionbox.cpp
	questionbox.h
	sortfunctions.cpp
	sortfunctions.h
	tableline.h
	termenergy.h
	vsolistelement.h
)

target_link_libraries(LineTable Qt5::Widgets)
