cmake_minimum_required(VERSION 3.1.0)

project(Calc)

set(CMAKE_AUTOMOC ON)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

add_library(Calc STATIC
    include/fit.h
	fit.cpp
)

target_include_directories(Calc PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(Calc PRIVATE
    Molecule
	Residualfit
	Tools
	Qt5::Widgets
)
