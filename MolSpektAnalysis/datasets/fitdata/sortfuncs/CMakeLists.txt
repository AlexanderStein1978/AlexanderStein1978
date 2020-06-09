cmake_minimum_required(VERSION 3.1.0)

project(FitDataSortFuncs)

set(CMAKE_AUTOMOC ON)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

add_library(FitDataSortFuncs STATIC
    sortbyltabandprogfunctor.cpp
	sortbyltabandprogfunctor.h
	sortforextractneworchangedwithoutsourcefunctor.cpp
	sortforextractneworchangedwithoutsourcefunctor.h
	sortfunctions.cpp
	sortfunctions.h
	tablelinesortfunctor.cpp
	tablelinesortfunctor.h
)

target_link_libraries(FitDataSortFuncs Qt5::Widgets)
