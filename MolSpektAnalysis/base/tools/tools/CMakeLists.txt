cmake_minimum_required(VERSION 3.1.0)

project(Tools)

set(CMAKE_AUTOMOC ON)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

add_library(Tools STATIC
	constants.h
    tools.cpp
	tools.h
	utils.cpp
	utils.h
)

target_link_libraries(Tools Qt5::Widgets)
