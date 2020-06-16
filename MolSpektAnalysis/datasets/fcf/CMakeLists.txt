cmake_minimum_required(VERSION 3.1.0)

project(FCF)

set(CMAKE_AUTOMOC ON)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

add_library(FCF STATIC
	include/fcfjdep.h
	include/fcftab.h
	include/fcftable.h
    fcfjdep.cpp
	fcftab.cpp
	fcftable.cpp
)

target_include_directories(FCF PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(FCF
	TableWindow
	Wavefuncs
	Qt5::Widgets
)
