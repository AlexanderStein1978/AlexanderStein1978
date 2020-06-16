cmake_minimum_required(VERSION 3.1.0)

project(TermView)

set(CMAKE_AUTOMOC ON)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

add_library(TermView STATIC
	include/termview.h
	termview.cpp
)

target_include_directories(TermView PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(TermView
	TableWindow
	Qt5::Widgets
)
