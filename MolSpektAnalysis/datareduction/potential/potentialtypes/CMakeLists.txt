cmake_minimum_required(VERSION 3.1.0)

project(PotentialTypes)

set(CMAKE_AUTOMOC ON)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

add_library(PotentialTypes STATIC
	include/AnaPot.h
	include/MLRPot.h
	include/MTTPot.h
	include/splinepot.h
	include/TangToenniesPot.h
    AnaPot.cpp
	MLRPot.cpp
	MTTPot.cpp
	splinepot.cpp
	TangToenniesPot.cpp
)

target_include_directories(PotentialTypes PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(PotentialTypes
	Potential
	Qt5::Widgets
)
