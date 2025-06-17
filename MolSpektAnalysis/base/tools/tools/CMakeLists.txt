cmake_minimum_required(VERSION 3.1.0)

project(Tools)

set(CMAKE_AUTOMOC ON)

find_package(Qt6 COMPONENTS Widgets REQUIRED)

add_library(Tools STATIC
	include/constants.h
	include/tools.h
	include/utils.h
    tools.cpp
	utils.cpp
)

target_include_directories(Tools PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(Tools
	Qt6::Widgets
)
