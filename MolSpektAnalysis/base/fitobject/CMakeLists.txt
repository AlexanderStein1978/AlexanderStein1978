cmake_minimum_required(VERSION 3.1.0)

project(FitObject)

set(CMAKE_AUTOMOC ON)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

add_library(FitObject STATIC
	include/fitobject.h
	include/gaussian.h
	include/lorentzian.h
    fitobject.cpp
	gaussian.cpp
	lorentzian.cpp
)

target_include_directories(FitObject PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(FitObject
	Tools
	Qt5::Widgets
)
