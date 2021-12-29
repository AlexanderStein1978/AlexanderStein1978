cmake_minimum_required(VERSION 3.1.0)

project(MdiChild)

set(CMAKE_AUTOMOC ON)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

add_library(MdiChild STATIC
	include/mdichild.h
    mdichild.cpp
)

target_include_directories(MdiChild PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(MdiChild PRIVATE
	MolSpektAnalysisMainWindow
	Qt5::Widgets
)
