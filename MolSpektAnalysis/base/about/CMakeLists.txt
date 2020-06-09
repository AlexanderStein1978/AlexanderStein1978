cmake_minimum_required(VERSION 3.1.0)

project(AboutMolSpektAnalysis)

set(CMAKE_AUTOMOC ON)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

add_library(AboutMolSpektAnalysis STATIC
    about.cpp
	about.h
)

target_link_libraries(AboutMolSpektAnalysis Qt5::Widgets)
