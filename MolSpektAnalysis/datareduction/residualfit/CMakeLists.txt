cmake_minimum_required(VERSION 3.1.0)

project(Residualfit)

set(CMAKE_AUTOMOC ON)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

add_library(Residualfit STATIC
	BandConst.h
    ConnectSpline.cpp
	ConnectSpline.h
	Js.h
	LocalPerturbation.cpp
	LocalPerturbation.h
	naturalspline.cpp
	naturalspline.h
	point.h
	ResidualFit.cpp
	ResidualFit.h
	Spline.cpp
	Spline.h
	SplinePoint.h
)

target_link_libraries(Residualfit Qt5::Widgets)
