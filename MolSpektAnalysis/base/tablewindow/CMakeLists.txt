cmake_minimum_required(VERSION 3.5.0)

project(TableWindow)

set(CMAKE_AUTOMOC ON)

find_package(Qt6 COMPONENTS Widgets REQUIRED)

add_library(TableWindow STATIC
	include/basedata.h
	include/tableviewwindow.h
	include/tableviewwindowcore.h
	include/tableviewwindowcoresortfunctor.h
	include/tablewidgetwindow.h
	include/tablewindow.h
	include/viewlist.h
	include/mtable.h
	include/tabsortfunctor.h
    mtable.cpp
    tableviewwindow.cpp
    tableviewwindowcore.cpp
    tablewidgetwindow.cpp
	tablewindow.cpp
)

target_include_directories(TableWindow PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(TableWindow PRIVATE
	MdiChild
	MolSpektAnalysisMainWindow
<<<<<<< HEAD
	Qt6::Widgets
=======
	Qt5::Widgets
>>>>>>> f91263c093dfe7215b3249af2e1113f12a7a6877
)
