cmake_minimum_required(VERSION 3.1.0)

project(MolSpektAnalysis)

set(CMAKE_AUTOMOC ON)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

add_executable(MolSpektAnalysis
    inputdialog.cpp
	inputdialog.h
	main.cpp
	MainWindow.cpp
	MainWindow.h
	showaction.h
)

target_link_libraries(MolSpektAnalysis Qt5::Widgets)
