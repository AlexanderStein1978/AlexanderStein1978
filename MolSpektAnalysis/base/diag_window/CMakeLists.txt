cmake_minimum_required(VERSION 3.1.0)

project(DiagWindow)

set(CMAKE_AUTOMOC ON)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

add_library(DiagWindow STATIC
	include/DiagWindow.h
	include/Picture.h
    DiagWindow.cpp
	Picture.cpp
)

target_include_directories(DiagWindow PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(DiagWindow
    MdiChild
	Tools
	Qt5::Widgets
)
