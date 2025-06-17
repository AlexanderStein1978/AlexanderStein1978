cmake_minimum_required(VERSION 3.1.0)

project(PotentialStructs)

set(CMAKE_AUTOMOC ON)

find_package(Qt6 COMPONENTS Widgets REQUIRED)

add_library(PotentialStructs STATIC
	include/elevel.h
	include/potentialdata.h
	include/wantedcoeffvalue.h
    elevel.cpp
	potentialdata.cpp
)

target_include_directories(PotentialStructs PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(PotentialStructs PRIVATE
	MonteCarloSim
	Tools
	Residualfit
	Qt6::Widgets
)
