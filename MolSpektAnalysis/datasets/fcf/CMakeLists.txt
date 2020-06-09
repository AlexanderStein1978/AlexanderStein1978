cmake_minimum_required(VERSION 3.1.0)

project(FCF)

set(CMAKE_AUTOMOC ON)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

add_library(FCF STATIC
    fcfjdep.cpp
	fcfjdep.h
	fcftab.cpp
	fcftab.h
	fcftable.cpp
	fcftable.h
)

target_link_libraries(FCF Qt5::Widgets)
