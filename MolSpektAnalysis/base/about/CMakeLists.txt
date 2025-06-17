cmake_minimum_required(VERSION 3.1.0)

project(AboutMolSpektAnalysis)

set(CMAKE_AUTOMOC ON)

find_package(Qt6 COMPONENTS Widgets REQUIRED)

add_library(AboutMolSpektAnalysis STATIC
	include/about.h
    about.cpp
)

target_include_directories(AboutMolSpektAnalysis PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(AboutMolSpektAnalysis
	Tools
	Qt6::Widgets
)
