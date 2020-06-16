cmake_minimum_required(VERSION 3.1.0)

project(DunhamCoefficients)

set(CMAKE_AUTOMOC ON)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

add_library(DunhamCoefficients STATIC
	include/duntable.h
    duntable.cpp
)

target_include_directories(DunhamCoefficients PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(DunhamCoefficients
	Calc
	LineTable
	Tools
	Qt5::Widgets
)
