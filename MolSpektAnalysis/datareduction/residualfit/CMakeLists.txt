cmake_minimum_required(VERSION 3.1.0)

project(Residualfit)

set(CMAKE_AUTOMOC ON)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

add_library(Residualfit STATIC
	include/BandConst.h
	include/Js.h
	include/LocalPerturbation.h
	include/naturalspline.h
	include/point.h
	include/ResidualFit.h
	include/SplinePoint.h
	include/Spline.h
    ConnectSpline.cpp
	ConnectSpline.h
	LocalPerturbation.cpp
	naturalspline.cpp
	ResidualFit.cpp
	Spline.cpp
)

target_include_directories(Residualfit PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(Residualfit PRIVATE
    Calc
	FitObject
	MolSpektAnalysisMainWindow
	Qt5::Widgets
)
