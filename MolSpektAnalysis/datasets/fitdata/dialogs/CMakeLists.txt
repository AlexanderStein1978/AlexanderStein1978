cmake_minimum_required(VERSION 3.1.0)

project(FitDataDialogs)

set(CMAKE_AUTOMOC ON)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

add_library(FitDataDialogs STATIC
	include/coupledpotfitimportdialog.h
	include/espotfitinputelstateassigndialog.h
	include/sourceoffsetdialog.h
	include/tsdialog.h
	include/vassigndialog.h
    coupledpotfitimportdialog.cpp
	espotfitinputelstateassigndialog.cpp
	sourceoffsetdialog.cpp
	tsdialog.cpp
	vassigndialog.cpp
)

target_include_directories(FitDataDialogs PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(FitDataDialogs
	MolSpektAnalysisMainWindow
	Qt5::Widgets
)
