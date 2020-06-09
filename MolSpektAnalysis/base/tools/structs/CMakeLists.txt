cmake_minimum_required(VERSION 3.1.0)

project(MolSpektAnalysisStructs)

set(CMAKE_AUTOMOC ON)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

add_library(MolSpektAnalysisStructs STATIC
    band.h
	fittetline.h
	fli.cpp
	fli.h
	line.h
)

target_link_libraries(MolSpektAnalysisStructs Qt5::Widgets)
