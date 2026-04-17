cmake_minimum_required(VERSION 3.5.0)

project(PotentialTypes)

set(CMAKE_AUTOMOC ON)

find_package(Qt6 COMPONENTS Widgets REQUIRED)

add_library(PotentialTypes STATIC
	include/AnaPot.h
	include/MLRPot.h
	include/MTTPot.h
	include/splinepot.h
	include/TangToenniesPot.h
    AnaPot.cpp
	MLRPot.cpp
	MTTPot.cpp
	splinepot.cpp
	TangToenniesPot.cpp
)

target_include_directories(PotentialTypes PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(PotentialTypes PRIVATE
	Calc
	LineTable
    MonteCarloSim
	Potential
	PotentialStructs
	Tools
<<<<<<< HEAD
	Qt6::Widgets
=======
	Qt5::Widgets
>>>>>>> f91263c093dfe7215b3249af2e1113f12a7a6877
)
