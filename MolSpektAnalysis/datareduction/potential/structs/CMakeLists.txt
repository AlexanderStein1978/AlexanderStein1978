cmake_minimum_required(VERSION 3.1.0)

project(PotentialStructs)

set(CMAKE_AUTOMOC ON)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

add_library(PotentialStructs STATIC
    elevel.cpp
	elevel.h
	potentialdata.cpp
	potentialdata.h
	wantedcoeffvalue.cpp
)

target_link_libraries(PotentialStructs Qt5::Widgets)
