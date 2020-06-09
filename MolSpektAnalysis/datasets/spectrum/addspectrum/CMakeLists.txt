cmake_minimum_required(VERSION 3.1.0)

project(LineTable)

set(CMAKE_AUTOMOC ON)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

add_library(LineTable STATIC
    adddialog.cpp
	adddialog.h
	addspektrum.cpp
	addspektrum.h
	cutdialog.cpp
	cutdialog.h
	doubleaddspektrum.cpp
	doubleaddspektrum.h
	markedpeak.h
)

target_link_libraries(LineTable Qt5::Widgets)
