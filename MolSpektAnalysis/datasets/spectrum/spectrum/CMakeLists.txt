cmake_minimum_required(VERSION 3.1.0)

project(Spektrum)

set(CMAKE_AUTOMOC ON)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

add_library(Spektrum STATIC
    exportlineprofiledialog.cpp
	exportlineprofiledialog.h
	feldialog.cpp
	feldialog.h
	selectlinetabledialog.cpp
	selectlinetabledialog.h
	spectlist.cpp
	spectlist.h
	Spektrum.cpp
	Spektrum.h
)

target_link_libraries(Spektrum Qt5::Widgets)
