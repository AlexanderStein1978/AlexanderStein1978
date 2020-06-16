cmake_minimum_required(VERSION 3.1.0)

project(TableWindow)

set(CMAKE_AUTOMOC ON)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

add_library(TableWindow STATIC
	include/tablewindow.h
	include/viewlist.h
	include/mtable.h
	include/tabsortfunctor.h
    mtable.cpp
	tablewindow.cpp
)

target_include_directories(TableWindow PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(TableWindow
	MdiChild
	Qt5::Widgets
)
