cmake_minimum_required(VERSION 3.1.0)

project(MolSpektAnalysisStructs)

set(CMAKE_AUTOMOC ON)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

add_library(MolSpektAnalysisStructs STATIC
    include/band.h
	include/fittetline.h
	include/fli.h
	include/line.h
	fli.cpp
)

target_include_directories(MolSpektAnalysisStructs PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(MolSpektAnalysisStructs Qt5::Widgets)
