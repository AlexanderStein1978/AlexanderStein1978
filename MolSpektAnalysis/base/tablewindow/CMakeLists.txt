cmake_minimum_required(VERSION 3.1.0)

project(TableWindow)

set(CMAKE_AUTOMOC ON)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

add_library(TableWindow STATIC
    mtable.cpp
	mtable.h
	tablewindow.cpp
	tablewindow.h
	tabsortfunctor.h
	viewlist.h
)

target_link_libraries(TableWindow Qt5::Widgets)
