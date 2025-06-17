cmake_minimum_required(VERSION 3.1.0)

project(FitDataSortFuncs)

set(CMAKE_AUTOMOC ON)

find_package(Qt6 COMPONENTS Widgets REQUIRED)

add_library(FitDataSortFuncs STATIC
	include/sortbyltabandprogfunctor.h
	include/sortforextractneworchangedwithoutsourcefunctor.h
	include/fitdatasortfunctions.h
	include/tablelinesortfunctor.h
    sortbyltabandprogfunctor.cpp
	sortforextractneworchangedwithoutsourcefunctor.cpp
	fitdatasortfunctions.cpp
	tablelinesortfunctor.cpp
)

target_include_directories(FitDataSortFuncs PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(FitDataSortFuncs
	MolSpektAnalysisMainWindow
	Qt6::Widgets
)
