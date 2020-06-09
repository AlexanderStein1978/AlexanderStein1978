cmake_minimum_required(VERSION 3.1.0)

project(FitDataDialogs)

set(CMAKE_AUTOMOC ON)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

add_library(FitDataDialogs STATIC
    coupledpotfitimportdialog.cpp
	coupledpotfitimportdialog.h
	espotfitinputelstateassigndialog.cpp
	espotfitinputelstateassigndialog.h
	sourceoffsetdialog.cpp
	sourceoffsetdialog.h
	tsdialog.cpp
	tsdialog.h
	vassigndialog.cpp
	vassigndialog.h
)

target_link_libraries(FitDataDialogs Qt5::Widgets)
