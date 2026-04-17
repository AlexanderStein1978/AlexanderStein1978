cmake_minimum_required(VERSION 3.5.0)

project(AddSpektrum)

set(CMAKE_AUTOMOC ON)

find_package(Qt6 COMPONENTS Widgets REQUIRED)

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
	Spektrum
<<<<<<< HEAD
	Qt6::Widgets
=======
	Qt5::Widgets
>>>>>>> f91263c093dfe7215b3249af2e1113f12a7a6877
)
