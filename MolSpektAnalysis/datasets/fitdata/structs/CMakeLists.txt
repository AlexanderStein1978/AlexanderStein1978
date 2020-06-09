cmake_minimum_required(VERSION 3.1.0)

project(FitDataStructs)

set(CMAKE_AUTOMOC ON)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

add_library(FitDataStructs STATIC
	assign_vs_comptrans.h
	levelcomb.h
	tlref.h
)

target_link_libraries(FitDataStructs Qt5::Widgets)
