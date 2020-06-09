cmake_minimum_required(VERSION 3.1.0)

project(PotentialTypes)

set(CMAKE_AUTOMOC ON)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

add_library(PotentialTypes STATIC
    AnaPot.cpp
	AnaPot.h
	MLRPot.cpp
	MLRPot.h
	MTTPot.cpp
	MTTPot.h
	splinepot.cpp
	splinepot.h
	TangToenniesPot.cpp
	TangToenniesPot.h
)

target_link_libraries(PotentialTypes Qt5::Widgets)
