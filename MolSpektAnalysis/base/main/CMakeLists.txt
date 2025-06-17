cmake_minimum_required(VERSION 3.1.0)

project(MolSpektAnalysis)

set(CMAKE_AUTOMOC ON)

find_package(Qt6 COMPONENTS Widgets REQUIRED)

add_executable(MolSpektAnalysis
	main.cpp
)

target_link_libraries(MolSpektAnalysis
	MolSpektAnalysisMainWindow
	Qt6::Widgets
)
