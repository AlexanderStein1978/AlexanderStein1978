cmake_minimum_required(VERSION 3.1.0)

project(FitDataStructs)

set(CMAKE_AUTOMOC ON)

find_package(Qt6 COMPONENTS Widgets REQUIRED)

add_library(FitDataStructs STATIC
	assign_vs_comptrans.h
	levelcomb.h
	tlref.h
)

target_include_directories(MolSpektAnalysisStructs PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(FitDataStructs Qt6::Widgets)
