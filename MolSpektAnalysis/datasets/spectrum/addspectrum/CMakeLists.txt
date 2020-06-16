cmake_minimum_required(VERSION 3.1.0)

project(AddSpektrum)

set(CMAKE_AUTOMOC ON)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

add_library(AddSpektrum STATIC
	include/adddialog.h
	include/cutdialog.h
    adddialog.cpp
	addspektrum.cpp
	addspektrum.h
	cutdialog.cpp
	doubleaddspektrum.cpp
	doubleaddspektrum.h
	markedpeak.h
)

target_include_directories(AddSpektrum PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(AddSpektrum
	DiagWindow
	Qt5::Widgets
)
